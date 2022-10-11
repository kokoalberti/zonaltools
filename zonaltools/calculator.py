
import types
import itertools
import math
import csv
import logging
import uuid
import os
import json

from collections import OrderedDict

import numpy as np
import numpy.ma as ma

from osgeo import gdal, ogr

logging.basicConfig(level=logging.DEBUG)

class NumpyEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.integer):
            return int(obj)
        if isinstance(obj, np.floating):
            return float(obj)
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        if isinstance(obj, np.bool_):
            return bool(obj)
        return super(NumpyEncoder, self).default(obj)

class Source(object):

    def __init__(self, source, name=None):
        self.source = source
        self.name = name

        self.ds = gdal.Open(self.source)

        if not self.ds:
            raise ValueError("Could not open datasource '{}' using GDAL.".format(self.source))

    def __repr__(self):
        return f"<Source name='{self.name}'>"


class Zone(object):

    def __init__(self, name=None, extras={}, identifier=None, wkt=None, feat=None, srs=None, geom_type=None, reference=None):
        self.name = name
        self.identifier = identifier
        self.reference = reference
        self.wkt = wkt
        self.feat = feat
        self.geom_type = geom_type
        self.srs = srs
        self.uuid = str(uuid.uuid4())
        self.cutline = "/vsimem/"+self.uuid+".fgb"
        self.cache = {}
        self.cutline_ds = None

        self.mask = None
        self.extras = extras
        self._lats = None
        self._lons = None
        self._area = None
        self._mem_paths = []

        self.results = {}
        self.custom_exports = OrderedDict()

    def __repr__(self):
        return f"<Zone name='{self.name}' identifier='{self.identifier}'>"

    @property
    def _export_directory(self):
        export_directory = "./exports/{}/".format(self.identifier)
        try: os.makedirs(export_directory)
        except: pass
        return export_directory

    def export(self, fieldnames=[]):
        export = OrderedDict()
        for fieldname in fieldnames:
            if fieldname in self.results:
                export[fieldname] = self.results[fieldname]
            if hasattr(self, fieldname):
                export[fieldname] = getattr(self, fieldname)
            if fieldname in self.extras:
                export[fieldname] = self.extras[fieldname]

        for item in self.custom_exports:
            export[item] = self.custom_exports[item]

        return export

    def export_property(self, name, value):
        logging.debug("exporting property: {}".format(name))
        self.custom_exports[name] = value


    def export_array(self, array, filename):
        logging.debug("Exporting array as GTiff...")
        logging.debug(array.dtype)

        if not filename.endswith(".tif"):
            raise ValueError("Only GeoTIFF exports supported at this time.")

        if array.shape != self.mask.shape:
            raise ValueError("You can only export arrays with the same shape as the mask.")


        gtiff = gdal.GetDriverByName("GTiff")
        gdal_datatype = None

        if array.dtype == np.byte:
            gdal_datatype = gdal.GDT_Byte

        if array.dtype == np.short:
            gdal_datatype = gdal.GDT_Int16

        if array.dtype == np.float32:
            gdal_datatype = gdal.GDT_Float32

        if array.dtype == np.float64:
            gdal_datatype = gdal.GDT_Float64

        logging.debug("Creating GTiff with GDAL datatype: {}".format(gdal_datatype))

        export_path = os.path.join(self._export_directory, filename)

        export_ds = gtiff.Create(export_path, xsize=array.shape[1], ysize=array.shape[0], bands=1, eType=gdal_datatype, options=['TILED=YES', 'COMPRESS=DEFLATE', 'ZLEVEL=2'])
        export_ds.SetGeoTransform(self.mask_geotransform)
        export_ds.SetProjection(self.mask_projection)
        export_ds.GetRasterBand(1).WriteArray(array)
        export_ds = None

        logging.debug("Done exporting!")


    def memfile(self, filename=None):
        if not filename:
            path = "/vsimem/{}.tif".format(str(uuid.uuid4()))
        else:
            path = "/vsimem/{}".format(filename)
        self._mem_paths.append(path)
        return path

    def unlink_all(self):
        for path in self._mem_paths:
            gdal.Unlink(path)

    def cleanup(self):
        self.unlink_all()
        del self.cache


    def create_mask(self):
        # TODO: Reading data here with Translate just in order to create the 
        #       mask. Can still be optimized.
        #
        # Get envelope of the feature that defines the zone
        logging.debug("Creating mask for zone {}".format(self))

        mem = ogr.GetDriverByName("Memory")
        poly_ds = mem.CreateDataSource("temp")
        poly_lyr = poly_ds.CreateLayer('polygons', None, self.geom_type)
        poly_lyr.CreateFeature(self.feat.Clone())


        co = ['TILED=YES', 'COMPRESS=DEFLATE', 'ZLEVEL=2']
        (minx, maxx, miny, maxy) = self.feat.GetGeometryRef().GetEnvelope()

        # Translate an envelope out of the reference source to do calculations on.
        template_ds = gdal.Translate(self.memfile(), self.reference.ds, projWin=(minx, maxy, maxx, miny), creationOptions=co)

        gtiff = gdal.GetDriverByName("GTiff")
        

        mask_path = self.memfile("mask.tif")

        mask_ds = gtiff.Create(mask_path, xsize=template_ds.RasterXSize, ysize=template_ds.RasterYSize, bands=1, eType=gdal.GDT_Byte, options=co)
        mask_ds.SetProjection(template_ds.GetProjection())
        mask_ds.SetGeoTransform(template_ds.GetGeoTransform())
        mask_ds.GetRasterBand(1).WriteArray(np.zeros((mask_ds.RasterYSize, mask_ds.RasterXSize), dtype=np.short))

        # Perform rasterization
        gdal.RasterizeLayer(mask_ds, [1], poly_lyr, burn_values=[1])

        self.mask_geotransform  = mask_ds.GetGeoTransform()
        self.mask_projection = mask_ds.GetProjection()
        self.mask = mask_ds.GetRasterBand(1).ReadAsArray().astype(np.byte)

        mask_ds.FlushCache()

        self.mask_ds = mask_ds
        self.mask_rows, self.mask_cols = self.mask.shape

        gt = self.mask_geotransform
        self.mask_xmin = min(gt[0], gt[0] + self.mask_cols * gt[1])
        self.mask_xmax = max(gt[0], gt[0] + self.mask_cols * gt[1])
        self.mask_ymin = min(gt[3], gt[3] + self.mask_rows * gt[5])
        self.mask_ymax = max(gt[3], gt[3] + self.mask_rows * gt[5])

        logging.debug("Pixel size of mask: x: {} y: {}".format(self.mask_geotransform[1], self.mask_geotransform[5]))
        logging.debug("Mask geotransform: {}".format(self.mask_geotransform))
        logging.debug("Done creating mask array of size {}x{} with {}/{} valid pixels.".format(self.mask.shape[0], self.mask.shape[1], self.mask.sum(), self.mask.size))

    def fetch_data(self, source):
        """
        Warp data from the source into a dataset matching the mask.
        """
        logging.debug("Fetching data for source '{}'...".format(source.name))
        outputBounds = (self.mask_xmin, self.mask_ymin, self.mask_xmax, self.mask_ymax)
        co = ['TILED=YES', 'COMPRESS=DEFLATE', 'ZLEVEL=2']

        data_ds = gdal.Warp(self.memfile(), 
                  source.source,
                  xRes=self.mask_geotransform[1],
                  yRes=self.mask_geotransform[5],
                  dstSRS=self.mask_ds.GetProjection(),
                  outputBounds=outputBounds,
                  creationOptions=co)

        data_arr = data_ds.GetRasterBand(1).ReadAsArray()
        data = ma.array(data_arr, mask=(self.mask==0))
        logging.debug("Done! Created masked array of size {}x{}.".format(data.shape[0], data.shape[1]))
        return data

    def lats(self):
        if self.mask is None:
            self.create_mask()

        if "lats" not in self.cache:
            lats_array = np.linspace(self.mask_ymax+(0.5*self.mask_geotransform[5]), self.mask_ymin-(0.5*self.mask_geotransform[5]), self.mask_rows)
            lats = ma.array(np.transpose([lats_array] * self.mask_cols), mask=(self.mask==0))
            self.cache["lats"] = lats

        return self.cache["lats"]

    def lons(self):
        if self.mask is None:
            self.create_mask()

        if "lons" not in self.cache:
            lons_array = np.linspace(self.mask_xmin+(0.5*self.mask_geotransform[1]), self.mask_xmax-(0.5*self.mask_geotransform[1]), self.mask_cols)
            lons = ma.array(np.tile(lons_array, (self.mask_rows, 1)), mask=(self.mask==0))
            self.cache["lons"] = lons

        return self.cache["lons"]

    def area(self):
        """
        Return a masked array containing area of the pixels in the mask.
        """
        if self.mask is None:
            self.create_mask()


        if "area" not in self.cache:
            pixel_size = ((-1*self.mask_geotransform[5])+self.mask_geotransform[1])/2
            lats_array = np.linspace(self.mask_ymax+(0.5*self.mask_geotransform[5]), self.mask_ymin-(0.5*self.mask_geotransform[5]), self.mask_rows)
            lats_area = [area_of_pixel(pixel_size, lat) for lat in lats_array]
            area = ma.array(np.transpose([lats_area] * self.mask_cols), mask=(self.mask==0))
            self.cache["area"] = area

        return self.cache["area"]

    def data_for_source(self, source):
        if self.mask is None:
            self.create_mask()

        if source.name not in self.cache:
            logging.debug("Fetching zone data (not cached) for '{}'".format(source.name))
            data = self.fetch_data(source)
            self.cache[source.name] = data
            return data
        else:
            logging.debug("Fetching zone data (cached) for '{}'".format(source.name))
            return self.cache[source.name]


class Calculator(object):
    """
    - auto fix projection differences
    - vectors as sources (either mask or extract field value and rasterize)
    - transparantly get/set custom properties to be exported
    - band access of data sources as "self.lossyear.b1" of self.lossyear[0]
    - performance review
    - generic export. if it's an array -> gtiff, if it's a dim0 -> csv
    - better management of export options 
    - run zones through a validate_zone() which can be overwritten by user for customizations
    - access to results per zone as an iterator
    - additional mask with weigting per pixel instead of 0/1 mask for more precision
    - provide a nodata mask representing nodata of the input file, and integrate this into the mask for the calculations to avoid NaNs in the final arrays
    - avoid as much as possile disk access, do everything in memory
    - implement the calculator as a with Calculator() as calc: calc.calculate() context to easily clean everything up when it's done.
    - consider implementing sources as class properties: sources = {'layer': 'layer.tif'} instead of (or together with add_source()
    - preloading remote datasource sources (/vsis3, /vsicurl etc file systems)
    - better error handling and debugging, esp. error "sre you sure you added a source with that name" when there is a syntax error
    - allow defining of other sources as the "master" copy to use for gridding
    - allow supplying your own grid to convert all the sources to
    - test suite (!!)
    - more reliable system for matching output data types generated to output tifs and csvs
    - output of geospatial properties???

    """

    # Define a list of public function names. When obtaining a list of user-
    # defined functions to calculate per zone we need this to easily identify
    # which are custom added, and which methods are part of the regular API.
    _public_function_names = ("calculate", "export", "add_zones", "add_source", "mask", "area", "lats", "lons", "export_array", "export_property", "before", "after")

    def __init__(self, debug=False):
        logging.debug("Initializing the calculator...")
        # Intialize some variables that we'll use to store zones, sources, and
        # the attributes we want to calculate.
        self._sources = {}
        self._zones = []
        self._attributes = list(self._get_custom_attributes())
        self._exportable_attributes = set()
        self._debug = debug
        self._custom_properties = dict()


        # Stores the zone that we're actively doing calculations on.
        self._active_zone = None

        self._reference_source_name = None

    def before(self):
        pass

    def after(self):
        pass

    def export_property(self, name, value):
        """
        Manually export a calculated property.

        Todo: solve this with getters/setters on own properties
        """
        self._exportable_attributes.add(name)
        self._active_zone.export_property(name, value)
        # logging.info("Expattr")
        # logging.info(self._exportable_attributes)



    def add_source(self, source, name=None):
        if name not in self._sources:
            self._sources[name] = Source(source, name=name)
            if not self._reference_source_name:
                logging.debug("Added source '{}' to the calculator. This is the reference source!!!".format(name))
                self._reference_source_name = name
            else:
                logging.debug("Added source '{}' to the calculator.".format(name))
        else:
            raise ValueError("You must supply a valid name for this source.")

    def add_zones(self, zones, name=None, layer=None, field=None, extra_fields=""):
        ds = ogr.Open(zones)

        extra_fields = extra_fields.split(",")

        if not ds:
            raise ValueError("Could not load zones from path: {}".format(zones))

        logging.debug("Loading zone features...")
        for layer in ds:
            logging.debug("Layer: {}".format(layer))
            layer_defn = layer.GetLayerDefn()
            layer_field_names = []
            for i in range(layer_defn.GetFieldCount()):
                layer_field_names.append(layer_defn.GetFieldDefn(i).GetName())

            if not field:
                field = layer_field_names[0]

            if field not in layer_field_names:
                raise ValueError("Could not find the field '{}' in the field names of this layer.".format(field))


            geom_type = layer.GetGeomType()
            srs = layer.GetSpatialRef()

            layer.ResetReading()
            for feature in layer:
                zone_id = feature.GetField(field)

                extras = {}
                for extra_field in extra_fields:
                    if extra_field in layer_field_names:
                        extras[extra_field] = feature.GetField(extra_field)
                        self._exportable_attributes.add(extra_field)

                logging.debug(" + Feature: {}".format(zone_id))
                geom = feature.geometry()

                z = Zone(name=name, identifier=zone_id, extras=extras, wkt=geom.ExportToWkt(), feat=feature.Clone(), geom_type=geom_type, srs=srs, reference=self._sources[self._reference_source_name])
                self._zones.append(z)

            # Only the first layer for now.
            break

    def export(self, filename=None, attributes=None, options={}, exclude=[]):

        if (filename and not filename.endswith('.csv')) and (filename and not filename.endswith('.json')):
            raise ValueError("Only CSV and JSON export currently supported.")

        logging.debug("Exporting results to file: {}".format(filename))
        logging.debug("Exportable attributes: {}".format(self._exportable_attributes))

        if filename == None:
            # Export as jsonable dict
            json_export = []
            fieldnames = ('identifier', 'name') + tuple(sorted(self._exportable_attributes))
            logging.debug("Exporting fieldnames: {}".format(fieldnames))
            for zone in self._zones:
                if zone.results or zone.custom_exports:
                    result = zone.export(fieldnames=fieldnames)
                    logging.debug("Result: {}".format(result))
                    json_export.append(result)

            # Return a encoded/decoded json to ensure that the serialization worked ok
            return json.loads(json.dumps(json_export, ensure_ascii=False, cls=NumpyEncoder).encode('utf-8'))


        if filename.endswith('.csv'):
            with open(filename, mode='w') as csv_file:
                fieldnames = ('identifier', 'name') + tuple(sorted(self._exportable_attributes)) + ('wkt', )
                writer = csv.DictWriter(csv_file, fieldnames=fieldnames, extrasaction='ignore', quoting=csv.QUOTE_ALL)
                writer.writeheader()
                for zone in self._zones:
                    if zone.results or zone.custom_exports:
                        logging.debug("Writing results for zone {}".format(zone))
                        writer.writerow(zone.export(fieldnames=fieldnames))

            return True

        if filename.endswith('.json'):
            json_export = []
            with open(filename, mode='wb') as json_file:
                fieldnames = ('identifier', 'name') + tuple(sorted(self._exportable_attributes))
                logging.debug("Exporting fieldnames: {}".format(fieldnames))
                for zone in self._zones:
                    if zone.results or zone.custom_exports:
                        result = zone.export(fieldnames=fieldnames)
                        logging.debug("Result: {}".format(result))
                        json_export.append(result)
                json_file.write(json.dumps(json_export, ensure_ascii=False, cls=NumpyEncoder).encode('utf-8'))

            return True

        logging.debug("Done!")

    def export_array(self, array, filename):
        """
        Exports an array as a GTiff file.
        """
        self._active_zone.export_array(array, filename)

    def _calculate_zone(self, zone):
        logging.debug("------------")
        logging.debug("Calculating data for a zone: {}".format(zone))
        logging.debug("Calculating the following attributes for this zone: {}".format(self._attributes))

        # Set the active zone.
        self._active_zone = zone

        # Run the before method
        self.before()

        for attribute in self._attributes:
            result = getattr(self, attribute)

            # If the attribute is actually callable, call it to fetch the 
            # results that it yielded. It must yield dicts, each of which 
            # updates the result dictionary.
            if callable(result):
                result_dict = {}
                for item in result():
                    if isinstance(item, dict):
                        result_dict.update(item)

                # Overwrite the original result variable
                result = result_dict

            dim = np.ndim(result)
            # Update the results dict for this attribute. Set to None in case
            # the dimension of the returned result is not 0, since those are 
            # most likely to be intermediary properties used to calculate other
            # values.
            if dim == 0:
                logging.debug("Dimension of result is 0!")
                # Anytime a value is returned with dimension 0 we add the name 
                # of that attribute to the "exportable_attributes" set so we 
                # know which attributes to export at the end of the run.
                self._exportable_attributes.add(attribute)
                self._active_zone.results.update({attribute: result})
            else:
                logging.debug("Dimension of result is not zero.")
                self._active_zone.results.update({attribute: None})

        logging.debug("Results for this zone:")
        logging.debug(self._active_zone.results)

        # Run the after method
        self.after()

        # Clear the cache on the zone we've just completed.
        self._active_zone.cleanup()
        logging.debug("------------")
        
        return True

    def calculate(self):
        """
        Calculate a property for a particular zone.
        """
        logging.debug("Loaded the following sources:")
        logging.debug(self._sources)
        logging.debug("Loaded the following zones:")
        logging.debug(self._zones)

        logging.debug("Starting calculations...")
        for zone in self._zones:
            if not zone.results:
                self._calculate_zone(zone)
        logging.debug("Done with calculations.")

    def _get_custom_attributes(self):
        """
        Returns a list of custom attributes that the user has defined on the 
        subclass.
        """
        for item in dir(self):
            if item.startswith("_"):
                continue
            if (item not in self._public_function_names):
                yield item

    @property
    def area(self):
        return self._active_zone.area()

    @property
    def lons(self):
        return self._active_zone.lons()

    @property
    def lats(self):
        return self._active_zone.lats()

    @property
    def mask(self):
        return self._active_zone.mask

    def __getattr__(self, name):
        """
        The __getattr__ method is used to let custom functions access data
        sources by name that have been defined using add_source().
        """
        if name in self._sources:
            source = self._sources[name]
            return self._active_zone.data_for_source(source)
        else:
            raise AttributeError(f"The attribute '{name}' is unknown. Did you add it as a source using 'add_source'?")


def zonal_property(f):
    @property
    def wrapper(self, *args, **kwargs):
        # if self.attribute is None:
        #     pass
        return f(self, *args, **kwargs)
    return wrapper

def area_of_pixel(pixel_size, center_lat):
    """
    From: https://gis.stackexchange.com/questions/127165/more-accurate-way-to-calculate-area-of-rasters

    Calculate m^2 area of a wgs84 square pixel.

    Adapted from: https://gis.stackexchange.com/a/127327/2397

    Parameters:
        pixel_size (float): length of side of pixel in degrees.
        center_lat (float): latitude of the center of the pixel. Note this
            value +/- half the `pixel-size` must not exceed 90/-90 degrees
            latitude or an invalid area will be calculated.

    Returns:
        Area of square pixel of side length `pixel_size` centered at
        `center_lat` in m^2.

    """
    a = 6378137  # meters
    b = 6356752.3142  # meters
    e = math.sqrt(1 - (b/a)**2)
    area_list = []
    for f in [center_lat+pixel_size/2, center_lat-pixel_size/2]:
        zm = 1 - e*math.sin(math.radians(f))
        zp = 1 + e*math.sin(math.radians(f))
        area_list.append(
            math.pi * b**2 * (
                math.log(zp/zm) / (2*e) +
                math.sin(math.radians(f)) / (zp*zm)))
    return pixel_size / 360. * (area_list[0] - area_list[1])
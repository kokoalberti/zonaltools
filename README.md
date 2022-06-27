# Zonaltools

Zonaltools is an experimental Python package to carry out custom zonal 
calculations on geospatial data sources.

Warning! This is an experimental package at the moment, many changes and 
features still need to be added and the API may change at any time.

# Example

Zonaltools works something like this:


```python
from zonaltools import Calculator, zonal_property

# Define the calculator
class MyCalculator(Calculator):

    # Anything decorated with zonal_property will be calculated for each zone
    # that is loaded into the calculator.
    @zonal_property
    def forest_area(self):
        """
        Calculate the forested area from the land use data source.
        """
        return ((self.landuse == 10) * self.area).sum()

# Initialize the calculator
calc = MyCalculator()

# Add the landuse.tif file as a source to the calculator using the 'add_source'
# method. The name specified in the 'name' kwarg is the attribute that the layer 
# can be accessed with in zonal calculations. So, in this case the landuse 
# raster data is accessible in the zonal calculation as 'self.landuse', which 
# will return a masked array of the data for that particular zone.
calc.add_source('landuse.tif', name='landuse')

# Add the zones that you want to run all the calculations. These can be vector
# files containing polygons. All the zonal properties defined in the calculator
# will be calculated for each zone in the source.
calc.add_zones('nationalparks.shp', name="nationalparks")

# Run the calculation
calc.calculate()

# Export the results. Will return the value of the 'forest_area' zonal property
# for each of the zones in nationalparks.shp.s
calc.export('output.json')

```

# Improvements

The following features and other improvements are still a work in progress:

* Tests
* Automatically reproject sources and zones in different CRSs.
* Ability to use vectors as sources as well. Calculator should rasterize them in the background.
* Band access of sources with more than one band
* Better code structure
* Installation via pip/PyPi
* More export options
* Better docs
* Weighted mask for pixels with fractional cover
* Nodata mask
* Preloading remote data sources
* Allow setting of options (e.g. preloading, grid selection) when initializing calculator
* Allow running calculator as context manager and have it clean up propertly when context is destroyed
* Allow defining of a custom grid to run the calculations on. Calculator should warp in the background.
* Allow setting a grid based on one of the sources
* Better error handling and debugging
* Better matching of output data types generated to the output file format (e.g. export an array as a tif, an int as a table)
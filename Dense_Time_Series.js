var visParams = {
  bands: ['B7', 'B5', 'B4'],
  min: 0,
  max: 5000,
  gamma: 0.6
};
Map.addLayer(Image1);
Map.addLayer(Image2);
//Get the image collection and filter it by path, row and date. Get an image collection of 1 year at a time
var collection1 = ee.ImageCollection('LANDSAT/LC08/C01/T1_SR').filter(ee.Filter.calendarRange(2016,2016,'year'))
.filter(ee.Filter.calendarRange(4,4,'month')).filter(ee.Filter.lt('CLOUD_COVER', 60))
.filter(ee.Filter.eq('WRS_PATH', 137))
.filter(ee.Filter.eq('WRS_ROW', 41)).first();
//Make an image collection from image 1 and 2. This is because we need to extract maximum NDVI of the two years
var NDVICollect = ee.ImageCollection.fromImages([Image1, Image2]);

//From the image collection made previously select Maximum NDVI and time corresponding to it
var Max_NDVI = NDVICollect.qualityMosaic('Max_NDVI');

//Add two Bands from Max_NDVI into the image of 2016, which will be used for logical filters
var Image1 = Image1.addBands(Max_NDVI.select('Max_NDVI').rename('NDVI2')).addBands(Max_NDVI.select('time2').rename('t2'))
.addBands(Max_NDVI.select('NBR2').rename('NBR3'));

//Create a logical conditions image for NDVI2 is greater than threshold NDVI, which is 0.5
var C1 = Image1.select('NDVI2').gt(0.5);


//Create a logical conditions image for NDVI2 being greater than NDVI1 by a threshold, which is 0.2
var C2 = ((Image1.select('NDVI2')).subtract(Image1.select('NDVI'))).gt(0.2);


//Create a logical condition image for NBR1 being greater than NBR2 by a threshold, which is 0.1
var C3 = ((Image2.select('NBR2').multiply(-1)).subtract(Image1.select('NBR').multiply(-1))).gt(0.1);


//Create a logical conditon image for the time of regeneration such that either the vegetation flourishing is earlier or after at least 150 days of burn
var C4 = (((Image1.select('t2')).subtract(Image1.select('time'))).gt(60)).or((Image1.select('time')).gt(Image1.select('t2')));

//Create a logical condition image for NBR1 being greater than NBR2 by a threshold, which is 0.1
var C5 = ((Image2.select('NBR').multiply(-1)).subtract(Image1.select('NBR').multiply(-1))).gt(0.1);

var C6 = (Image1.select('NBR')).eq(Image1.select('NBR2'));
Map.addLayer(C6, {}, 'C6');

//Multiply all the conditions
var conditions = C1.multiply(C2).multiply(C4).multiply(C3);
Map.addLayer(collection1, visParams);
Map.addLayer(Image1.select('classification').gt(0.5).multiply(conditions),
{min:0, max: 1, palette:['white', 'red'] });

var Fire = (Image1.select('classification').multiply(conditions)).gt(0.75);
var growingconditions = (Image1.select('classification').gt(0.35)).multiply(C1);//.multiply(C3);
Map.addLayer(Fire);
//Map.addLayer(Firms);
var mask = Fire.gt(0);
Map.addLayer(mask);
var boxcar = ee.Kernel.square({radius: 1, units: 'pixels'});
var pixcount = (Fire.updateMask(mask)).connectedPixelCount({maxSize: 100});
var mainfire = pixcount.gt(11);
var mainfire = mainfire.unmask();

var opened = mainfire.focal_max({kernel:boxcar});
var Fire1 = opened.multiply(growingconditions);
Map.addLayer(Fire1, {}, 'Fire1');
 for (var i= 1; i <= 200; i = i + 1){
  var Intermediate = Fire1.focal_max({kernel:boxcar});
  var Fire2 = Intermediate.multiply(growingconditions);
  var Fire1 = Fire2;
 }
Export.image.toDrive({
  image:Fire2,
  description:'Fire_Npl11',
  scale:30,
  //Assign region as the geometry of the first image of the collection
  region:Image1.geometry(),
  maxPixels: 1e11,
  //Assign the crs of the image as it appears on the original collection
  crs: 'EPSG:32645'
}); 


Map.addLayer(Fire2, {min:0, max: 1, palette:['white', 'red'] }, 'Fire2');
Map.addLayer(table);
Map.addLayer(C3);
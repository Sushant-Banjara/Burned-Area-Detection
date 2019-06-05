//***The script takes in Landsat 8 image collection of each year and applies random forest classification to output burned area probability of each pixel***//
//***Input 1 (Train_img): The image corresponding to April 2016, from which the training pixel is extracted***//
//***Input 2 (Fire): The polygon corresponnding to burned training pixels in Train_img***//
//***Input 3 (No_Fire): The polygon corresponding to unburned training pixels in Train_img***//
//***Input 4 (LC): Landcover map of Nepal for 2010 used to mask different landcover features, if required***//
//***Output: An image containing maximum burned area probability attained by a pixel in a given year and corresponding NDVI, NBR, and Time***//

//Get the name of the input bands and output bands
var inBands = ee.List(['B2','B3','B4','B5','B6','B7', 'IC']);
var outBands = ee.List(['blue','green','red','nir','swir1','swir2', 'IC']);

//Set visual paramaters for later visualization of images
var visParams = {
  bands: ['swir2', 'nir', 'green'],
  min: 0,
  max: 0.5,
  gamma: 0.6
};

//Get the image collection and filter it by path, row, date, and percentage cloud cover. Get the image collection of one year at a time
var collection1 = ee.ImageCollection('LANDSAT/LC08/C01/T1_SR').filter(ee.Filter.calendarRange(2015,2015,'year'))
.filter(ee.Filter.calendarRange(1,12,'month')).filter(ee.Filter.lt('CLOUD_COVER', 40))
.filter(ee.Filter.eq('WRS_PATH', 144))
.filter(ee.Filter.eq('WRS_ROW', 40));


//Function for masking out unnecessary and/or confounding pixels. The masking function uses bit shifting
function masks(image) {
  // Bits 3 and 5 are cloud shadow and cloud, respectively.
  var cloudShadowBitMask = (1 << 3);
  var cloudsBitMask = (1 << 5);
  //Bit 2 is mask for water and 4 is mask for snow
  var snowBitMask = (1 <<4);
  var waterBitMask = (1 << 2);
  // Get the pixel QualityAssurance band
  var qa = image.select('pixel_qa');
  // Use bitwise And operator to mask out pixels
  var mask = qa.bitwiseAnd(cloudShadowBitMask).eq(0)
                 .and(qa.bitwiseAnd(cloudsBitMask).eq(0)).and(qa.bitwiseAnd(waterBitMask).eq(0)).and(qa.bitwiseAnd(snowBitMask).eq(0));
      
  return image.updateMask(mask);
}

//Map the masking function to the entire collection so that unnecessary pixels in each image of the collection are removed
var collection = collection1.map(masks);
//////////////////////////////////////////////////////////////////Image Preprocessing Ends Here///////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//Get the Digital Elevation Model
var DEM = ee.Image("USGS/SRTMGL1_003");

// Define a boxcar or low-pass kernel
var boxcar = ee.Kernel.square({
radius: 1, units: 'pixels'});

//Smoothen the DEM using low-pass filter
var DEMs = DEM.convolve(boxcar);

//Create a constant value to convert degrees into radians
var rad = 0.01745;

//Compute the slope and aspect of each pixel and convert to radians
var slope = ee.Terrain.slope(DEMs).multiply(rad);
var aspect = ee.Terrain.aspect(DEMs).multiply(rad);

//Compute sine and cosine of the slope
var s1 = slope.cos();
var s2 = slope.sin();

//Function to compute the solar illumination condition of each pixel and add as a band to original image
//The variables names are the same as the variable naes used in the original paper
function IlluminationCondition(image){
  //Get the solar zenith angle of an image and compute its sine and cosine
  var z = ee.Image.constant(ee.Number(image.get('SOLAR_ZENITH_ANGLE'))).multiply(rad);
  var z1 = z.cos();
  var z2 = z.sin();
  
  //Get the solar azimuth angle of an image and compute its sine and cosine
  var a = ee.Image.constant(ee.Number(image.get('SOLAR_AZIMUTH_ANGLE'))).multiply(rad);
  
  var solar = s1.multiply(z1);
  var asp1 = s2.multiply(z2);
  var asp2 = (a.subtract(aspect)).cos();
  var asp = asp1.multiply(asp2);
  
  //Compute Illumination condition using variables computed previously
  var IC = solar.add(asp);
  
  return image.addBands(IC.rename('IC'));
}

//Map the IlluminationCondition function to the entire image collection
var collection = collection.map(IlluminationCondition);


//Select only the required bands and change default names to outBands for easy identification
var collection = collection.select(inBands, outBands);

//Create the list of the band that will be topographically corrected
var bandList = ['blue','green','red','nir','swir1','swir2'];

//Function to select each band from image and correct the band using statistical empirical method of topographic correction
var coll = function (img){
  
  //Time of an image is required for further analysis so don't let it get lost during image processing
  var st= img.get('system:time_start');
 
 //Function to calculate mean and linear regression for each band
 function mean(band){
  var x = img.select(band).reduceRegion({
        reducer: ee.Reducer.mean(),
        geometry: ee.Geometry(img.geometry()),
        scale: 30,
        maxPixels: 1e11
      });
  var y = img.select('IC', band).reduceRegion({
        reducer: ee.Reducer.linearFit(),
        geometry: ee.Geometry(img.geometry()),
        scale: 30,
        maxPixels: 1e11
      });
    
   var m = ee.Number(y.get('scale'));
   var b = ee.Number(y.get('offset'));
   var Ilavg = ee.Number(x.get(band));
   var corrected = img.expression(
        "(image - (b + m*ic) + mean)*0.0001", {
        'image': img.select(band),
        'ic': img.select('IC'),
        'b': b,
        'm': m,
        'mean': Ilavg}
        );
        return corrected.double();
}
  //Map the function to entire bandlist and return as an image
  return  ee.Image(bandList.map(mean)).set('system:time_start', st);
  
  
};

//Map the function of topographic correction to entire image collection
var corrected = collection.map(coll);
///////////////////////////////////////////////Topographic Correction Ends Here//////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Write the function to create the necessary spectral indices and add as band in the image
function SpectralIndices(image){
var SWIR1 = image.select('swir1');
var SWIR2= image.select('swir2');
var NIR = image.select('nir');
var Red = image.select('red');
var NDVI = image.normalizedDifference(['nir', 'red']).rename ('NDVI').double();
var NBR = (image.normalizedDifference(['nir', 'swir2'])).multiply(-1).rename('NBR').double();
var NBR2 = image.normalizedDifference(['swir1', 'swir2']).rename('NBR2').double(); 
var NDMI = image.normalizedDifference(['nir', 'swir1']).rename('NDMI').double();
var MIRBI = ((SWIR2.multiply(10)).subtract(SWIR1.multiply(9.8))).add(2);
var SAVI = ((NIR.subtract(Red)).multiply(1.5)).divide((NIR.add(Red)).add(0.5));
var BAI = (((NIR.subtract(0.06)).pow(2)).add((Red.subtract(0.1)).pow(2))).pow(-1);
var N = (((NIR.pow(2)).subtract(Red.pow(2))).multiply(2).add(NIR.multiply(1.5)).add(Red.multiply(0.5))).divide(NIR.add(Red).add(0.5)); 
var GEMI = (N.subtract((N.pow(2)).multiply(0.25)).subtract(Red.subtract(0.125))).divide((Red.multiply(-1)).add(1));

return image.addBands(NDVI.rename('NDVI')).addBands(NBR.rename('NBR')).addBands(NBR2.rename('NBR2')).addBands(NDMI.rename('NDMI'))
.addBands(MIRBI.rename('MIRBI')).addBands(SAVI.rename('SAVI')).addBands(BAI.rename('BAI')).addBands(GEMI.rename('GEMI'));
}


//Map Spectral Indices to entire image collection
var corrected = corrected.map(SpectralIndices);

//Selected and display one image from corrected image collection to check everything is okay
var Imagelist = corrected.toList(corrected.size());
var kimg= ee.Image(Imagelist.get(5));
Map.addLayer(kimg, visParams, 'kimg');
print(ee.Date(kimg.get('system:time_start')));
print(kimg);

//Merge the traing regions as Region of Interest
var ROI = Fire.merge(No_Fire);

//Create a random forest classifier
var classifier = ee.Classifier.randomForest({
	numberOfTrees: 100,
});

//Set the output of the classifier to probability
var classifier = classifier.setOutputMode('PROBABILITY');

//Create a list of bands used in classification
var bands = ['blue','green','red','nir','swir1','swir2', 'NDVI', 'NBR', 'NBR2', 'NDMI', 'MIRBI', 'SAVI', 'BAI', 'GEMI']; 

//Apply SpectralIndices function to Training Image
var Trainer = SpectralIndices(Train_img);

//Define the sample region from which the training data is extracted
var training = Trainer.sampleRegions({
	collection: ROI,
	properties: ['Class'],
	scale: 30
});

//Train the Model
var train = classifier.train(training, 'Class', bands);

//Write a function for classifying rest of the image from collection
function classification (kmg){
 var classified = kmg.classify(train);
 var classified1 = classified.double();
 var time_start = kmg.metadata('system:time_start');
  var base_band = kmg.select('red');
  
  //Change from UNIX time to days and add as band
  var time_band = base_band.add(
      time_start.divide(1000 * 60 * 60 * 24)
    ).rename('time');
 return classified1.addBands(time_band).addBands(kmg.select('NDVI')).addBands(kmg.select('NBR')); 
}

//Map the classification function to entire collection
var fire_prob = corrected.map(classification);
/////////////////////////////////////////////////////////Burned Area Classification Ends Here//////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//Create an image that has the Maximum NDVI within the collection. This will later be added on the image to be exported
var Max_NDVI = fire_prob.qualityMosaic('NDVI');
print(Max_NDVI);

//Create an image of Maximum NBR (actually it's minimum NBR but we multiplied NBR with negative 1 for ease)
var Max_NBR = fire_prob.qualityMosaic('NBR');


//From the entire collection select pixels with highest probability of being fire
var highest_prob = fire_prob.qualityMosaic('classification');


//Add the maximum NDVI, corresponding time, and maximum NBR to the image that has to be exported
var highest_prob = highest_prob.addBands(Max_NDVI.select('NDVI').rename('Max_NDVI')).addBands(Max_NDVI.select('time').rename('time2'))
.addBands(Max_NBR.select('NBR').rename('NBR2'));

//Print the property of the image
print(highest_prob);

//Visualize the highest_prob map to check everything works okay
Map.addLayer(highest_prob);

//Export image to Asset
Export.image.toAsset({
  image:highest_prob,
  description:'2015_144_40',
  scale:30,
  //Assign region as the geometry of the first image of the collection
  region:collection1.first().geometry(),
  maxPixels: 1e11,
  //Assign the crs of the image as it appears on the original collection
  crs: 'EPSG:32644'
});
////////////////////////////////////////////////The First Part of Algorithm Ends Here//////////////////////////////////////
///////////For Next Part Similar Image of Previous Year is produced for Time Series Analysis and Rejection Tests///////////
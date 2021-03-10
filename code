//Purpose:
//calculate standard deviation of Bands 1-7 minus sd of Bands 4-7
//then mask water, shadow, snow, cloud, get maximum value composite
//export final product

var exportGeometry = ee.Geometry.Rectangle(-63.2792, -9.6807, -60.6644, -11.6832);

//Uncomment one desired location
//Map.setCenter(-0.529897, 6.32609, 16); //Ghana cloud area

//Map.setCenter(-0.5737, 6.332, 12); //Ghana

Map.setCenter(-62.414, -10.554, 8); //Brazil

var bounds = Map.getBounds(true)
// Make a date filter to get images in this date range.
var dateFilter = ee.Filter.date('1999-01-01', '2013-12-31');

//specify bad frames in this format
var badFrames = ee.List(['LE07_193056_20130207'])

// Load SR for band data
//var SR = ee.ImageCollection('LANDSAT/LE7_SR')
var SR = ee.ImageCollection('LANDSAT/LE07/C01/T1_SR')
    .filter(dateFilter)
    .filter(ee.Filter.inList('system:index',badFrames).not());
    
var L7_2000 = SR.filterDate('1999-01-01', '2003-12-31');
var L7_2009 = SR.filterDate('2009-05-01', '2009-09-30');
// Creates three extra bands within an Image:
//  - stdev_B1-B7: The standard deviation of Bands 1, 2, 3, 4, 5, & 7 (no B6)
//  - stdev_B4-B7: The standard deviation of Bands 4, 5, & 7 (no B6)
//  - stdev_diff: The difference between stdev_B1-B7 & stdev_B4-B7
var stdev_diff = function(image) {
  var stdev = image
    .addBands(
      image.select("B1", "B2", "B3", "B4", "B5", "B7").reduce(ee.Reducer.stdDev())
        .rename("stdev_B1-B7")).round()
    .addBands(
      image.select("B4", "B5", "B7").reduce(ee.Reducer.stdDev())
       .rename("stdev_B4-B7").round());
  return stdev
    .addBands(stdev.select("stdev_B1-B7").subtract(stdev.select("stdev_B4-B7"))
    .rename("stdev_diff"));
}
// This function converts QA info from binary - need to update
var qualityBit = function(image) {
  var quality = image.select("pixel_qa")
  var pixelQA = quality.bitwiseAnd(0x2).rightShift(1).rename('clear');
  return image.addBands(pixelQA);
};

//Function to mask water, shadow, snow, cloud - only clear pixels remain, then apply 5x5 focal min
var maskingmin = function(image){
  var cfmask = image.select('clear').eq(1);
  var kernel = ee.Kernel.square({
      radius: 75,
      units: 'meters',
    });
  var mincfmask = cfmask.reduceNeighborhood(
      ee.Reducer.min(), kernel);
  return image.updateMask(mincfmask);
 // return image.updateMask(cfmask.and(B1).and(B2).and(B3));
}


//This creates the maximum value composite after cloud masking
var SVVI_2000_maxfm = L7_2000.map(qualityBit).map(stdev_diff).map(maskingmin).select('stdev_diff').max()
var SVVI_2009_maxfm = L7_2009.map(qualityBit).map(stdev_diff).map(maskingmin).select('stdev_diff').max()

//This function applies the 3x3 focal SD
var focalSD = function(img) {
  var kernel = ee.Kernel.square({
      radius: 45,
      units: 'meters',
    });
  return img.reduceNeighborhood(
      ee.Reducer.stdDev(), kernel);
};
var f3x3_2000 = L7_2000.map(qualityBit).map(stdev_diff).map(maskingmin).select('stdev_diff').map(focalSD).max();
var f3x3_2009 = L7_2009.map(qualityBit).map(stdev_diff).map(maskingmin).select('stdev_diff').map(focalSD).max();

// apply the 5x5 focal mean
var f5x5_2000 = f3x3_2000.focal_mean(75,'square','meters'); //Display varies with resolution
var f5x5_2009 = f3x3_2009.focal_mean(75,'square','meters'); //Display varies with resolution
//Adjust thresholds here
var f5x5threshold = 40
var SVVIthreshold = 150
//applies a threshold to 5x5 for values greater than 40
var f5x5_2000_threshold = f5x5_2000.mask(f5x5_2000.lte(f5x5threshold));
var f5x5_2009_threshold = f5x5_2009.mask(f5x5_2009.lte(f5x5threshold));
//unmask so that values of 0 will be used properly in map algebra
var f5x5_2000_thresholdUM = f5x5_2000.lte(f5x5threshold).unmask();
var f5x5_2009_thresholdUM = f5x5_2009.lte(f5x5threshold).unmask();
//applies threshold to values greater than 150 in original SVVI image
var SVVI_2000_threshold = SVVI_2000_maxfm.mask(SVVI_2000_maxfm.lte(SVVIthreshold));
var SVVI_2009_threshold = SVVI_2009_maxfm.mask(SVVI_2009_maxfm.lte(SVVIthreshold));
//unmask for map algebra
var SVVI_2000_thresholdUM = SVVI_2000_maxfm.lte(SVVIthreshold).unmask();
var SVVI_2009_thresholdUM = SVVI_2009_maxfm.lte(SVVIthreshold).unmask();


var forest_2000 = f5x5_2000_thresholdUM.and(SVVI_2000_thresholdUM);
var forest_2009 = f5x5_2009_thresholdUM.and(SVVI_2009_thresholdUM);

//Subract forest layers
var forestDiff = forest_2000.subtract(forest_2009)//loss = 1, gain = -1, no change = 0


//mask layers
var forest_2000only = forest_2000.updateMask(forest_2000.eq(1))
var forest_2009only = forest_2009.updateMask(forest_2009.eq(1))
//update mask so that areas of no change are not shown
var forestDiffonly = forestDiff.updateMask(forestDiff.eq(0).not())

//Mask and convert binary

var SVVI_2009_thresholdMask = SVVI_2009_threshold.gt(SVVIthreshold)

//These visulization parameters change the bands displayed and the range shown   
var vizParams = {  bands: ['stdev_diff'], min: -300,  max: 500,};

//add layers to map

Map.addLayer(SVVI_2009_maxfm, vizParams, "2009 GEE SVVI continuous"); 
Map.addLayer(SVVI_2009_thresholdMask, vizParams, "2009 Thresholded SVVI"); 
Map.addLayer(exportGeometry,{},"area to export")

Export.image.toDrive({
  image: SVVI_2009_maxfm,
  description: 'SVVI_2009_maxfm',
  scale: 30,
  maxPixels:1e12 ,
  region: exportGeometry
})

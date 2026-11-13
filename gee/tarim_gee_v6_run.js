// ============================================================
//  Lower Tarim River - NDVI / LST series 2000-2023        (script v6)
//  Run at https://code.earthengine.google.com
//
//  What changed with respect to v5
//  ------------------------------
//  v5 buffered a centreline that had been traced on the NDVI field, i.e. the
//  green ribbon, not the river.  v6 buffers the actual mapped Tarim channel
//  (OpenStreetMap waterway=river, stitched from Daxihaizi to Taitema), so the
//  corridor is now "within 1 km of the river channel" in the literal sense.
//
//  Channel: 43 vertices, 375.1 km of centreline, 1 km buffer on both sides.
//  Everything else (products, season, scaling, seed, statistics) is unchanged
//  from v5, so the new export can be compared with the old one value by value.
//
//  Run the four exports ONE AT A TIME (Export > Run, wait for "Completed").
//  With all of them queued at once the earlier run failed.
// ============================================================

// ---------- 0. settings ---------------------------------------------------
var YEARS    = ee.List.sequence(2000, 2023);
var MON_FROM = 4, MON_TO = 10;      // growing season, April ... October
var BUFFER_M = 1000;                // corridor half width, metres (0-1 km)
var FITBOX   = ee.Geometry.Rectangle([86.2, 38.9, 90.4, 41.5]);  // unchanged
var SEED     = 42;

var NDVI_ID = 'MODIS/061/MOD13Q1';  // 250 m, 16 day
var LST_ID  = 'MODIS/061/MOD11A2';  // 1 km,   8 day

// ---------- 1. river channel --------------------------------------------
// OSM waterway=river centreline, Daxihaizi Reservoir -> Taitema Lake.
// Quote these coordinates in Section 2.1 and archive them with the code.
var channel = ee.Geometry.LineString([
  [87.3795,40.6474], [87.4323,40.6255], [87.4828,40.6064],
  [87.5508,40.5879], [87.6019,40.5655], [87.6480,40.5557],
  [87.6843,40.5269], [87.7357,40.5428], [87.7842,40.5213],
  [87.8235,40.4923], [87.8628,40.4613], [87.9056,40.4417],
  [87.9525,40.4321], [88.0010,40.4441], [88.0465,40.4321],
  [88.1043,40.4345], [88.1457,40.4068], [88.1757,40.3716],
  [88.2177,40.3448], [88.2560,40.3146], [88.2583,40.2752],
  [88.2696,40.2347], [88.2576,40.1984], [88.2916,40.1652],
  [88.3413,40.1591], [88.3759,40.1328], [88.3738,40.0887],
  [88.3892,40.0442], [88.3803,40.0025], [88.3782,39.9570],
  [88.3813,39.9178], [88.4169,39.8878], [88.3944,39.8462],
  [88.4144,39.8022], [88.4097,39.7571], [88.4378,39.7175],
  [88.4377,39.6753], [88.4491,39.6242], [88.4377,39.5848],
  [88.4298,39.5483], [88.4143,39.5127], [88.3859,39.4754],
  [88.3685,39.4709]
]);

var corridor = channel.buffer(BUFFER_M);

// ---------- 2. composites -------------------------------------------------
function ndviYear(y) {
  var start = ee.Date.fromYMD(y, MON_FROM, 1);
  var end   = ee.Date.fromYMD(y, MON_TO, 1).advance(1, 'month');
  return ee.ImageCollection(NDVI_ID).select('NDVI')
    .filterDate(start, end).mean().multiply(0.0001).rename('NDVI');
}

function lstYear(y) {
  var start = ee.Date.fromYMD(y, MON_FROM, 1);
  var end   = ee.Date.fromYMD(y, MON_TO, 1).advance(1, 'month');
  var k = ee.ImageCollection(LST_ID).select('LST_Day_1km')
    .filterDate(start, end).mean().multiply(0.02).rename('LST_K');
  return k.addBands(k.subtract(273.15).rename('LST'));      // deg C
}

// NDVI keeps its native 250 m; LST is resampled onto the same grid.
function stackYear(y) {
  var nd = ndviYear(y);
  var ls = lstYear(y).resample('bilinear')
    .reproject({ crs: nd.select('NDVI').projection(), scale: 250 });
  return nd.addBands(ls);
}

// ---------- 3. corridor statistics ---------------------------------------
function statsYear(y) {
  var a = ndviYear(ee.Number(y)).reduceRegion({
    reducer: ee.Reducer.mean().combine({
      reducer2: ee.Reducer.percentile([50, 90]), sharedInputs: true }),
    geometry: corridor, scale: 250, maxPixels: 1e13 });
  var b = lstYear(ee.Number(y)).reduceRegion({
    reducer: ee.Reducer.mean(), geometry: corridor, scale: 1000, maxPixels: 1e13 });
  return ee.Feature(null, { year: y,
    NDVI_mean: a.get('NDVI_mean'), NDVI_p50: a.get('NDVI_p50'),
    NDVI_p90: a.get('NDVI_p90'), LST_C: b.get('LST') });
}

// ---------- 4. SELF-CHECK - read the console before exporting ------------
// Expect (1) length ~375 km, area ~750 km2;  (2) NDVI roughly 0.05-0.10 at
// the start of the record and 0.08-0.13 at the end.  A value below 0.04 means
// the buffer is sitting on the bare channel bed or on desert.
print('channel length [km]', ee.Number(channel.length()).divide(1000));
print('corridor area  [km2]', corridor.area().divide(1e6));
print('corridor WKT', corridor.toGeoJSONString());
print('QA - six years (read this before you export)',
      ee.FeatureCollection(ee.List([2000, 2005, 2010, 2015, 2020, 2023])
        .map(function (y) { return statsYear(y); })));

// ---------- 5. map -------------------------------------------------------
Map.centerObject(channel, 8);
Map.addLayer(ee.ImageCollection(NDVI_ID).select('NDVI')
  .filterDate('2000-01-01', '2024-01-01')
  .filter(ee.Filter.calendarRange(MON_FROM, MON_TO, 'month'))
  .mean().multiply(0.0001), { min: 0.0, max: 0.35,
  palette: ['white', 'yellow', 'green', 'darkgreen'] }, 'mean NDVI 2000-2023');
Map.addLayer(channel,  { color: 'yellow' }, 'OSM channel', true);
Map.addLayer(corridor, { color: 'red' },    'corridor = 1 km buffer', false);
Map.addLayer(FITBOX,   { color: 'blue' },   'TVDI fitting box', false);

// ---------- 6a. EXPORT 1 - annual series --------------------------------
var annual = ee.FeatureCollection(YEARS.map(function (y) { return statsYear(y); }));
print('annual series (also exported)', annual);
Export.table.toDrive({ collection: annual,
  description: 'tarim_annual_ndvi_lst_2000_2023_series', fileFormat: 'CSV' });

// ---------- 6b. EXPORT 2 and 3 - point tables for TVDI -------------------
// Comment these out for the first run, uncomment after export 1 completes.
var ptsCorridor = ee.FeatureCollection.randomPoints(corridor, 3000, SEED);
var ptsRegion   = ee.FeatureCollection.randomPoints(FITBOX, 3000, SEED);

function sampleYear(y, pts, zone) {
  var stack = stackYear(ee.Number(y))
    .addBands(ee.Image.constant(y).toFloat().rename('year'));
  return stack.reduceRegions({ collection: pts, reducer: ee.Reducer.mean(),
                               scale: 250, tileScale: 4 })
    .map(function (f) {
      var c = f.geometry().coordinates();
      return f.set({ zone: zone, lon: c.get(0), lat: c.get(1) });
    });
}

Export.table.toDrive({ collection: ee.FeatureCollection(
    YEARS.map(function (y) { return sampleYear(y, ptsCorridor, 'corridor'); })).flatten(),
  description: 'tarim_points_corridor_2000_2023', fileFormat: 'CSV',
  selectors: ['year', 'zone', 'NDVI', 'LST', 'LST_K', 'lon', 'lat'] });

Export.table.toDrive({ collection: ee.FeatureCollection(
    YEARS.map(function (y) { return sampleYear(y, ptsRegion, 'region'); })).flatten(),
  description: 'tarim_points_region_2000_2023', fileFormat: 'CSV',
  selectors: ['year', 'zone', 'NDVI', 'LST', 'LST_K', 'lon', 'lat'] });

// ---------- 7. OPTIONAL - distance bands (near-water effect) -------------
// Run this only after the three exports above are done.  It gives the mean
// NDVI of the 0-0.5, 0.5-1, 1-2 and 2-3 km bands, which is the evidence for
// the "benefit decays with distance from the channel" statement.
/*
function bandMean(y, inner, outer) {
  var g = channel.buffer(outer).difference(channel.buffer(inner));
  var m = ndviYear(ee.Number(y)).reduceRegion({ reducer: ee.Reducer.mean(),
    geometry: g, scale: 250, maxPixels: 1e13 });
  return ee.Feature(null, { year: y, inner: inner, outer: outer,
                            NDVI: m.get('NDVI') });
}
var bands = [];
[[0, 500], [500, 1000], [1000, 2000], [2000, 3000]].forEach(function (b) {
  YEARS.getInfo().forEach(function (y) { bands.push(bandMean(y, b[0], b[1])); });
});
Export.table.toDrive({ collection: ee.FeatureCollection(bands),
  description: 'tarim_ndvi_distance_bands', fileFormat: 'CSV' });
*/

// ---------- 8. methods wording -------------------------------------------
// "The riparian corridor was defined once, independently of the year, as the
// 1 km buffer on each side of the mapped Tarim channel between the Daxihaizi
// Reservoir (87.3795 E, 40.6474 N) and Taitema Lake (88.3685 E, 39.4709 N);
// the centreline is 375 km long and the buffer covers 750 km2.  NDVI is the
// growing-season (April-October) mean of MOD13Q1 at 250 m, LST that of
// MOD11A2 at 1 km, both 2000-2023."

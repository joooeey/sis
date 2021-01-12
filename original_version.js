/*
MIT License

Copyright (c) 2019 Lukas Schreiber at the University of British Columbia

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

print('script started');

// =======================================
// options
// =======================================

// smelly global variable
var max_ew_slope = 0; // below this east-west slope value, both ascending and descending are used. Null means no filtering.
// constants
var OBSERVATION_BOUNDS = 1826; // days before today that can be selected with sliders
var NOW = new Date(); // time of script start
var MIN_OBS = 3; // min number of datapoints required at right tail of step function

// =======================================
// ascending/descending separation based on terrain
// =======================================

// get aspect and slope bands
DEM = ee.Terrain.products(DEM);
  
// get relief image with Eastern sun
var aspect = DEM.select('aspect').multiply(Math.PI/180);
var slope = DEM.select('slope').multiply(Math.PI/180);
var ew_slope_img = slope.tan().multiply(aspect.sin()).atan().multiply(180/Math.PI);

// function to filter steep faces facing the satellite
function without_steep (img) {
  // asc: >0 (truthy), desc: 0 (falsy)
  var ascending = ee.String(img.get('orbitProperties_pass')).compareTo('DESCENDING');
  return ee.Algorithms.If(ascending,
                          img.updateMask(ew_slope_img.gt(-max_ew_slope)),
                          img.updateMask(ew_slope_img.lte(max_ew_slope))
  // for max_ew_slope == 0, desc orbit will be used for flat areas (i. e. where relief == 0)
                         );
}

// =======================================
// get season model
// =======================================
// code block adapted from Karyn Tabor:
// https://code.earthengine.google.com/26b5ce2fe826f4fe1d6b030cc229b07e

var EPOCH = ee.Date('2010-01-01');
function add_independents(image) {
  var freq_x_time = image.date().difference(EPOCH, 'year').multiply(2*Math.PI);
  var independents = ee.Image.cat(1, freq_x_time.sin(), freq_x_time.cos()).rename(['center', 'sin', 'cos']);
  return image.addBands(independents).float().select(['center', 'sin', 'cos', 'sigma0']);
}

var season_model;
function get_season_model(collection) {
  var with_independents = collection.map(add_independents); 
  var regression = with_independents.reduce(ee.Reducer.linearRegression(3), 4)
                                    .select('coefficients')
                                    .arrayProject([0])
                                    .arrayFlatten([['center', 'sin', 'cos']])
                                    .float();

  season_model = regression;
}

// function to remove seasonal component at each pixel.
var remove_season = function(image) {

  var freq_x_time = ee.Date(image.get('system:time_start')).difference(EPOCH, 'year').multiply(2*Math.PI);

  return ee.Image(image.subtract(season_model.select('sin').multiply(freq_x_time.sin()))
                  .subtract(season_model.select('cos').multiply(freq_x_time.cos()))
                  .float()
                  .copyProperties(image, ['system:time_start']));

};

// =======================================
// change visualization functions
// =======================================

// single step function (least squares fit)
function fit_step(full_coll, change_coll, min_obs) {

  var combi_reducer = ee.Reducer.mean().combine(
                      ee.Reducer.variance(), '', true).combine(
                      ee.Reducer.count(), '', true);

  var lsq_sums = change_coll.map( function(image) {
    image = ee.Image(image);
    full_coll = ee.ImageCollection(full_coll);
    change_coll = ee.ImageCollection(change_coll);
    var step_date = image.get('system:time_start');
    var left_tail = full_coll.filterMetadata('system:time_start', 'less_than', step_date);
    var left_stats = left_tail.reduce(combi_reducer);
    var right_tail = change_coll.filterMetadata('system:time_start', 'not_less_than', step_date);
    var right_stats = right_tail.reduce(combi_reducer);
    // add up two tails
    var neg_squares = left_stats.expression("b('sigma0_variance')*b('sigma0_count')")
                      .add(right_stats.expression("b('sigma0_variance')*b('sigma0_count')"))
                      .multiply(-1)
                      .where(right_stats.select('sigma0_count').lte(min_obs), -1e9);
    return ee.Image.cat(image.metadata('system:time_start'), neg_squares,
                        left_stats.select('sigma0_mean'), right_stats.select('sigma0_mean'))
              .rename(['toc', 'good_fit', 'sigma0_old', 'sigma0_new']);
  });

  var mosaic = ee.ImageCollection(lsq_sums).qualityMosaic('good_fit');
  mosaic = mosaic.updateMask(mosaic.select('good_fit').neq(-1e9)); // remove pixels were constraints are too hard.

  // R squared:
  var blue = mosaic.expression("1 + b('good_fit')/(stats.sigma0_variance*stats.sigma0_count)",
                               {stats: full_coll.reduce(ee.Reducer.variance().combine(
                                                          ee.Reducer.count(), '', true))
                               })
    .rename('R_squared');
  var green = mosaic.select('sigma0_new');
  var red = mosaic.select('sigma0_old');
  var full_count = full_coll.count().rename('images_in_time_series');
  var change_count = change_coll.count().rename('images_in_change_series');
  return ee.Image.cat(red, green, blue, full_count, change_count);
}

// trend line (least squares fit)
function fit_line(change_coll, start, end, weights) {
  change_coll = ee.ImageCollection(change_coll);
  start = ee.Date(start).millis();
  end = ee.Date(end).millis();
  ///*
  var set_weights;
  if (weights == 'no_weights') set_weights = function (image) { return image; };
  else if (weights == 'linear') set_weights = function (image) {
    return image.updateMask(ee.Number(image.get('system:time_start')).subtract(start).divide(end.subtract(start)));
  };
  else throw 'Bad argument supplied for weights!';
  change_coll = change_coll
    .map(function(image) {
      return set_weights(image)
        .addBands(image.metadata('system:time_start'));
    } );
  var trend_img = change_coll
    .select(['system:time_start', 'sigma0'])
    .reduce(ee.Reducer.linearFit());
  if (weights != 'no_weights') trend_img = trend_img.updateMask(1);

  // change offset (=intercept) from offset in 1970 to offset 'a_while_ago' and adjust units 
  var output_img = ee.Image.cat(
    trend_img.select('offset')
      .add(trend_img.select('scale').multiply(start)), // move intercept
    trend_img.select('offset')
      .add(trend_img.select('scale').multiply(end)), // move intercept
    // r_squared (actually nothing)
    ee.Image.constant(0),
    // count
    change_coll.select('sigma0').count()
  ).rename(['before', 'after', 'r_2', 'count']).float();
  return output_img;
}

// This function highlights backscatter maxima in color
// hue: date of maximum, saturation: sigma0 ratio max/median, value: median [decibels]
function find_spike(change_coll, start, end) {
  
  // main computations
  var spike_images = change_coll.map(function(image) {
    return image.addBands(image.metadata('system:time_start'));
  });
  var spike_mosaic = spike_images.qualityMosaic('sigma0');
  var median = spike_images.select('sigma0').median().rename('baseline');
  var diff = spike_mosaic.select('sigma0').subtract(median).rename('difference');
  var ratio = diff.expression('10.0**(b()/10.0)'); // =sigma0_max/sigma0_median
  
  // compute the dates for the legend
  start = ee.Date(start).millis().getInfo();
  var range = ee.Date(end).millis().getInfo() - start;
  var dates = [{bg:'red', c:'white', d:start},
               {bg:'yellow', c:'black', d: start + range/5},
               {bg:'00FF00', c:'white', d: start + range*2/5},
               {bg:'cyan', c:'black', d: start + range*3/5},
               {bg:'blue', c:'white', d: start + range*4/5},
               {bg:'magenta', c:'white', d: start + range}
              ];
  var legends = [ui.Label('spike date (topmost spike layer): ')];
  for (var i=0; i<dates.length; i++) {
    legends.push(ui.Label(ee.Date(dates[i].d).format('yyyy-MM-dd').getInfo(),
                          {color: dates[i].c, backgroundColor: dates[i].bg}
                         ));
  }
  var colorbar = ui.Panel(legends, ui.Panel.Layout.flow('horizontal'), {position: 'bottom-center'});
  Map.add(colorbar);
  
  // create hsv summary image
  var rgb = ee.Image.cat(spike_mosaic.select('system:time_start').subtract(start).divide(range*6/5), // clipping to hues from red to magenta, i.e. rygcbm
                         ratio.subtract(1).divide(9), // clipping to ratios from 1 to 10.
                         median.add(20).divide(20)
                        ).hsvToRgb();
  return rgb;
}

// =======================================
// export basic functions
// =======================================

exports.without_steep = without_steep;
exports.get_season_model = get_season_model;
exports.remove_season = remove_season;
exports.fit_line = fit_line;
exports.fit_step = fit_step;
exports.find_spike = find_spike;

// =======================================
// custom date slider for UI
// =======================================

// factory function
function single_slider (slider_default) {

  // define UI elements
  var mini_panel = ui.Panel({
    widgets: [
        ui.Slider({
          min: -OBSERVATION_BOUNDS,
          max: 1,
          value: slider_default,
          step: 1,
          style: {
            width: '200px',
            fontSize: 0, // hide label
            margin: '0px -50px 0px 8px',
          }
        }),
        ui.Textbox({style: {
          width: '85px',
          margin: '0px'
        } })
      ],
    layout: ui.Panel.Layout.flow('horizontal'),
    style: {margin: '0px'},
    }
  );
  
  // helper functions with convenient names
  var listen_textbox = true;
  var listen_slider = true;
  mini_panel.slider = mini_panel.widgets().get(0);
  mini_panel.box = mini_panel.widgets().get(1);
  mini_panel.box.parent = mini_panel;
  mini_panel.value = mini_panel.slider.getValue; // a function
  mini_panel.getDate = function() {
    mini_panel.box.style().set('color', 'black');
    if (listen_slider) {
      var date = new Date(NOW);
      date.setUTCDate(date.getUTCDate() + mini_panel.value());
      var pad = function(int) { return ('0' + int).slice(-2); };
      var date_string = date.getUTCFullYear() + '-' + pad(date.getUTCMonth()+1) + '-' + pad(date.getUTCDate());
      listen_textbox = false;
      mini_panel.box.setValue(date_string);
      listen_textbox = true;
      return date;
    }
  };
  
  // link textbox and slider
  mini_panel.slider.onSlide(mini_panel.getDate);
  mini_panel.box.onChange( function(text, box) {
    if (listen_textbox) {
      var newdate = new Date(text);
      var newval = Math.ceil((newdate - NOW)/(1000*3600*24));
      if (newval === 0 || newval && newval > -OBSERVATION_BOUNDS && newval < 2) {
        listen_slider = false;
        box.parent.slider.setValue(newval);
        listen_slider = true;
      }
      else
        box.style().set('color', 'red');
    }
  });
  
  // sync textbox and slider
  mini_panel.getDate();
  
  // return combined slider with textbox
  return mini_panel;
}

// =======================================
// UI
// =======================================

// help box
var help_text = ui.Label({
  value: 'This application can summarize a Sentinel-1 timeseries in one RGB image. The options on the right side control how. Do not forget to zoom in! The "inspect pixels" toggle on the left side allows you to plot the time series at any pixel.',
  style: {position: 'bottom-center', width: '370px', whiteSpace: 'pre-wrap'},
});
var help_cross = ui.Button({
  label: 'X',
  style: {position: 'top-right'},
});
var help_panel = ui.Panel({
  layout: ui.Panel.Layout.absolute(),
  widgets: [help_cross, help_text],
  style: {width: '400px', height: '200px'},
});
Map.add(help_panel);
help_cross.onClick( function() {help_panel.style().set('shown', false); });
function show_help_panel(text) {
  help_panel.style().set('shown', true);
  help_text.setValue(text);
}

// toggle to remove harmonic (seasonal) signal
var season_toggle = ui.Checkbox({
  label: 'remove seasonal signal',
  style: {stretch: 'both'},
});
var season_help = ui.Button({
  label: '?',
  onClick: function() {
    show_help_panel('At each pixel a sine-function is fit to the timme series of pixel values. Subsequently this seasonal component is subtracted across the time series. Note that this step only makes sense if there is at least a year between start of data and start of monitoring period. Even then this step tends to introduce noise.');
  }
});
var season_panel = ui.Panel({
  widgets: [season_toggle, season_help],
  layout: ui.Panel.Layout.flow('horizontal'),
});

// dropdown to select polarization
var pol_label = ui.Label('polarization');
var pol_dropdown = ui.Select({
  items: ['VH', 'VV'],
  value: 'VV',
  style: {stretch: 'horizontal'},
});
var pol_help = ui.Button({
  label: '?',
  onClick: function() {
    show_help_panel('VH: Transmitted radio waves are vertically polarized and horizontally polarized waves are received.\n\nVV: Both transmitted and received radio waves are vertically polarized.');
  }
});
var pol_panel = ui.Panel({
  widgets: [pol_dropdown, pol_help],
  layout: ui.Panel.Layout.flow('horizontal'),
  style: {stretch: 'horizontal'},
});

// textbox to select layer name
var name_label = ui.Label('layer name:');
var name_box = ui.Textbox({
  value: 'step composite image',
  style: {stretch: 'horizontal'},
});

// dropdown to select between ascending and descending polarization or combinations thereof
var pass_label = ui.Label('satellite pass direction:');
var pass_dropdown = ui.Select({
  items: ['ascending', 'descending', 'best', 'filter steep terrain (>20°)', 'combine'],
  value: 'combine',
  style: {stretch: 'horizontal'},
});
var pass_help = ui.Button({
  label: '?',
  onClick: function() {
    switch (pass_dropdown.getValue()) {
      case 'ascending':
        show_help_panel('Only images from ascending satellite passes are used to build the composite image.');
        break;
      case 'descending':
        show_help_panel('Only images from descending satellite passes are usedto build the composite image.');
        break;
      case 'best':
        show_help_panel('For each pixel, only images from satellites passing on the uphill side are used.');
        break;
      case 'filter steep terrain (>20°)':
        show_help_panel('For pixels in steep terrain (>20°), only images from satellites passing on the uphill side are used.');
        break;
      case 'combine':
        show_help_panel('The full time series with both ascending and descending node images is used to build the composite image.');
        break;
    }
  }
});
var pass_panel = ui.Panel({
  widgets: [pass_dropdown, pass_help],
  layout: ui.Panel.Layout.flow('horizontal'),
  style: {stretch: 'horizontal'},
});

// dropdown to select time series fitting function
var fit_label = ui.Label('function to fit:');
var fit_dropdown = ui.Select({
  items: ['trend line', 'step function', 'spike function'],
  value: 'step function',
  style: {stretch: 'horizontal'},
  onChange: function(new_value) {
    name_box.setValue(new_value);
    d1_label.style().set('shown', (new_value == 'step function'));
    d1_slide.style().set('shown', (new_value == 'step function'));
  },
});
var fit_help = ui.Button({
  label: '?',
  onClick: function() {
    switch (fit_dropdown.getValue()) {
      case 'trend line':
        show_help_panel('At each pixel a linear fit of backscatter values is computed. The composite shows the start (red) and end (green) values of this trend line.');
        break;
      case 'step function':
        show_help_panel('At each pixel a step function with a single discontinuity is fit to the backscatter values. The composite shows the left (red) and right (green) values of the step function, as well as R squared in blue.');
        break;
      case 'spike function':
        show_help_panel('At each pixel the maximum and average of the time series is computed. The resulting composite shows the following. value: average; saturation: max/average; hue: Date of maximum.');
        break;
    }
  }
});
var fit_panel = ui.Panel({
  widgets: [fit_dropdown, fit_help],
  layout: ui.Panel.Layout.flow('horizontal'),
  style: {stretch: 'horizontal'},
});

// sliders to select extent of monitoring period
var d1_label = ui.Label('start of data:');
var d1_slide = single_slider(-180);
var d2_label = ui.Label('start of monitoring period:');
var d2_slide = single_slider(-90);
var d3_label = ui.Label('end of monitoring period:');
var d3_slide = single_slider(0);

// button to calculate
var apply_button = ui.Button({
  label: 'load time series composite image',
  style: {stretch: 'horizontal'},
});

// panel with all the options
var side_panel = ui.Panel({
  widgets: [season_panel,
    pol_label, pol_panel, // not implemented
    pass_label, pass_panel,
    fit_label, fit_panel, name_label, name_box,
    d1_label, d1_slide, d2_label, d2_slide, d3_label, d3_slide,
    apply_button],
});
ui.root.add(side_panel);

// function for loading the time series composite images
var show = function() {
  print('loading ' + fit_dropdown.getValue() + ' composite image');
  
  var start_date = d1_slide.getDate();
  var steady_date = d2_slide.getDate();
  var end_date = d3_slide.getDate(); 
  
  print(ee.Date(start_date), ee.Date(steady_date), ee.Date(end_date));
  
  // preprocessing
  var Sen1filtered = Sen1.filterDate(start_date, end_date)
                         .filter(ee.Filter.listContains('transmitterReceiverPolarisation', pol_dropdown.getValue()))
                         .select(pol_dropdown.getValue()) // not implemented
                         .map(function (img) {return img.rename('sigma0')});
  
  switch (pass_dropdown.getValue()) {
    case 'ascending':
      Sen1filtered = Sen1filtered.filter(ee.Filter.eq('orbitProperties_pass', 'ASCENDING'));
      break;
    case 'descending':
      Sen1filtered = Sen1filtered.filter(ee.Filter.eq('orbitProperties_pass', 'DESCENDING'));
      break;
    case 'best':
      max_ew_slope = 0;
      Sen1filtered = Sen1filtered.map(without_steep);
      break;
    case 'filter steep terrain (>20°)':
      max_ew_slope = 20;
      Sen1filtered = Sen1filtered.map(without_steep);
      break;
    case 'combine':
      break;
  }

  var Sen1filtered_steady = Sen1filtered.filterDate(start_date, steady_date);

  if (season_toggle.getValue()) {
    // get harmonic season model and remove it from images
    get_season_model(Sen1filtered_steady);
    Sen1filtered = Sen1filtered.map(remove_season);
    Sen1filtered_steady = Sen1filtered.filterDate(start_date, steady_date);
  }
  var Sen1timeframe = Sen1filtered.filterDate(steady_date, end_date);

  var viz_rgb;
  switch (pol_dropdown.getValue()) {
    case 'VV':
      viz_rgb = {
        min: [-20.0, -20.0, 0.0],
        max: [0.0, 0.0, 1.0],
      };
      break;
    case 'VH':
      viz_rgb = {
        min: [-30.0, -30.0, 0.0],
        max: [0.0, 0.0, 1.0],
      };
      break;
  }
  var layer_name = name_box.getValue();
  switch (fit_dropdown.getValue()) {
    case 'trend line':
      var trend_composite = fit_line(
        Sen1timeframe, steady_date, end_date, 'no_weights');
      Map.addLayer(trend_composite, viz_rgb, layer_name);
      break;
  case 'step function':
      var step_composite = fit_step(
        Sen1filtered, Sen1timeframe, MIN_OBS);
      Map.addLayer(step_composite, viz_rgb, layer_name);
    break;
  case 'spike function':
      var widgets = Map.widgets();
      if (widgets.length() > 2) Map.remove(widgets.get(2)); // remove previous colorbar
      var spike_img = find_spike(
        Sen1timeframe, steady_date, end_date);
      Map.addLayer(spike_img, {}, layer_name);
    break;
  }
};
apply_button.onClick(show);

// =======================================
// Export Map Button
// =======================================

var export_map_button = ui.Button({
  label: 'export a layer',
  style: {position: 'top-right'},
  onClick: function (button) {
    var layers = Map.layers();
    var layer_objects = [];
    layers.forEach(function (layer, index) {
      layer_objects.push({label: layer.getName(), value: layer});
    });
    layer_objects.push({label: 'export all', value: layers});
    Map.add(ui.Select({
      items: layer_objects,
      style: {position: 'top-right'},
      onChange: function(layer, dropdown) {
        if (layer instanceof ui.Map.Layer) {
          Export.image.toDrive({
            image: layer.getEeObject().float().set('visualization', layer.getVisParams()),
            description: layer.getName().replace(/[^a-z0-9]/gi, '_'),
            folder: 'new_change_composites',
            scale: 10, // meters
            region: Map.getBounds(true),
            maxPixels: 1e11 // increase pixel limit to 1e11 pixels, allowing bigger files
          });
        }
        else if (layer instanceof ui.data.ActiveList) {
          layers.forEach(function (layer, index) {
            Export.image.toDrive({
              image: layer.getEeObject().float().set('visualization', layer.getVisParams()),
              description: layer.getName().replace(/[^a-z0-9]/gi, '_'),
              folder: 'new_change_composites' + new Date().toISOString(),
              scale: 10, // meters
              region: Map.getBounds(true),
              maxPixels: 1e11 // increase pixel limit to 1e11 pixels, allowing bigger files
            });
          });
        }
        else throw 'Invalid selection for export!';
        Map.remove(dropdown);
      }
    }));
  }
});

Map.add(export_map_button);

// =======================================
// Show time series on click
// =======================================

// separate ascending and descending columns
function asc_des_columns (img) {
  var band_names = img.bandNames().remove('angle');
  var new_names = band_names.map(function(name) {
    return ee.String(name).cat('_').cat(ee.String(img.get('orbitProperties_pass')).toLowerCase());
  });
  return img.select(band_names, new_names);
}
var Sen1_asc_des = Sen1.map(asc_des_columns);

// chart options
var vAxis_dB = {title: 'SAR backscatter coefficient sigma-0 [dB]',
                minValue: -30.0,
                maxValue: 0.0,
                gridlines: {count: 7}};
var options = {
  title: 'SAR backscatter time series',
  vAxis: vAxis_dB,
  hAxis: {title: 'acquisition start dates (UTC)'},
  lineWidth: 0,
  pointSize: 2,
  colors: ['#9999ff', '#ff9999', '#0000aa', '#aa0000'],
};

// checkbox: if checked clicking the map will produce a time series chart
function change_cursor(bool, box) {
  if(bool) Map.style().set({cursor: 'crosshair'});
  else Map.style().set({cursor: 'hand'});
}
var box_inspect = ui.Checkbox({
  label: 'show pixel time series',
  value: false,
  onChange: change_cursor,
  style: {backgroundColor: 'green', color: 'white', padding: '10px',
          stretch: 'horizontal', position: 'top-left'}
});
Map.add(box_inspect);

// function to compute time series chart
function show_chart (coords, map) {
  if(box_inspect.getValue()) { 
    // get point from coordinates
    var coord_array = Object.keys(coords).map(function (key) { return coords[key]; });
    var point = ee.Geometry.Point(coord_array);
    // create chart
    options.title = 'SAR backscatter time series at lon: ' + String(coords.lon) + ' lat: ' + coords.lat;
    var chart = ui.Chart.image.series(Sen1_asc_des, point)
                .setOptions(options);
    chart.style().set({width: '900px', height: '280px', position: 'bottom-left',
                       padding: '0px', margin: '0px',
    });
    // figure out if click is on bottom or top half of the map
    var center_coords = Map.getCenter();
    var y_center = center_coords.coordinates().get(1);
    var position = (coords.lat > y_center.getInfo()) ? 'bottom-center' : 'top-center';
    // show chart
    var chart_panel = ui.Panel({
      widgets: [chart,
                ui.Button({
                  label: 'X',
                  style: {position: 'top-right',
                          padding: '0px', margin: '0px',
                  },
                  onClick: function () { Map.remove(chart_panel); },
                }),
      ],
      layout: ui.Panel.Layout.absolute(),
      style: {width: '1000px', height: '300px', position: position,
              padding: '0px', margin: '0px',
      },
    });
    Map.add(chart_panel);
  }
}
Map.onClick(show_chart);

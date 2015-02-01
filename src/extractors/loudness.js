
function Loudness(opts) {
  if (!(this instanceof Loudness)) return new Loudness(opts);

  // featureInfo for loudness
  this.info = {
    "type": "multipleArrays",
    "arrayNames": {
      "1": "total",
      "2": "specific"
    }
  };

  this.NUM_BARK_BANDS = opts.NUM_BARK_BANDS || 24;
  this.normalisedSpectrum = opts.normalisedSpectrum;
  this.sampleRate = opts.sampleRate;

  this.specific = new Float32Array(this.NUM_BARK_BANDS);
  this.process = this.process.bind(this); // optimize later
  this.barkScale = opts.barkScale;
  this.bbLimits  = this.computeBarkBandLimits(this.barkScale, this.normalisedSpectrum.length, this.NUM_BARK_BANDS);
}

Loudness.prototype.computeBarkBandLimits = function(barkScale, nSpectrumLength, num_bark_bands) {
  barkScale       = barkScale       || this.barkScale;
  nSpectrumLength = nSpectrumLength || this.normalisedSpectrum.length;
  num_bark_bands  = num_bark_bands  || this.NUM_BARK_BANDS;
  
  var currentBandEnd = barkScale[nSpectrumLength-1]/num_bark_bands;
  var currentBand = 1;

  var bbLimits = new Int32Array(num_bark_bands+1);
  bbLimits[0] = 0;

  for(var i = 0; i<nSpectrumLength; i++){
   while(barkScale[i] > currentBandEnd) {
     bbLimits[currentBand++] = i;
     currentBandEnd = currentBand*barkScale[nSpectrumLength-1]/num_bark_bands;
   }
  }

  bbLimits[num_bark_bands] = nSpectrumLength-1;

  return bbLimits;
};

Loudness.prototype.computeSpecificLoudness = function(bbLimits, specific, nSpectrum, num_bark_bands) {
  bbLimits        = bbLimits       || this.bbLimits;
  specific        = specific       || this.specific;
  nSpectrum       = nSpectrum      || this.normalisedSpectrum;
  num_bark_bands  = num_bark_bands || this.NUM_BARK_BANDS;

  // console.log(bbLimits, specific, nSpectrum, num_bark_bands);

  for (var i = 0; i < num_bark_bands; i++){
    var sum = 0;

    for (var j = bbLimits[i] ; j < bbLimits[i+1] ; j++) {
      sum += nSpectrum[j];
    }
    console.log(sum)
    specific[i] = Math.pow(sum, 0.23);
  }

  return specific;
};

// belongs in utils
Loudness.prototype.sumArray = function(array) {
  var sum = 0;

  for (var i = 0; i < array.length; i++) {
    sum += array[i];
  }

  return sum;
};

Loudness.prototype.process = function() {
  var barkScale = this.barkScale;
  var nSpectrum = this.normalisedSpectrum;
  var nSpectrumLength = nSpectrum.length;
  var NUM_BARK_BANDS = this.NUM_BARK_BANDS;
  var specific = this.specific;
  var bbLimits  = this.bbLimits;

  //process
  var spec  = this.computeSpecificLoudness(bbLimits, specific, nSpectrum, NUM_BARK_BANDS);
  // console.log(spec, this.specific);
  var total = this.sumArray(spec);

  return {
    specific: spec,
    total: total
  };
};

module.exports = Loudness;
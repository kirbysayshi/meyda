
function Loudness(opts) {
  if (!(this instanceof Loudness)) return new Loudness(opts);

  this.NUM_BARK_BANDS = NUM_BARK_BANDS = opts.NUM_BARK_BANDS || 24;
  this.normalisedSpectrum = opts.normalisedSpectrum;
  this.sampleRate = opts.sampleRate;

  this.specific = new Float32Array(NUM_BARK_BANDS);
  this.process = this.process.bind(this);
  this.barkScale = opts.barkScale;
}

Loudness.prototype.createBarkBandLimits = function() {
  var nSpectrumLength = this.normalisedSpectrum.length;
  var NUM_BARK_BANDS = this.NUM_BARK_BANDS;
  var barkScale = this.barkScale;
  var bbLimits = new Int32Array(NUM_BARK_BANDS+1);
  
  bbLimits[0] = 0;
  var currentBandEnd = barkScale[nSpectrumLength-1]/NUM_BARK_BANDS;
  var currentBand = 1;
  for(var i = 0; i<nSpectrumLength; i++){
   while(barkScale[i] > currentBandEnd) {
     bbLimits[currentBand++] = i;
     currentBandEnd = currentBand*barkScale[nSpectrumLength-1]/NUM_BARK_BANDS;
   }
  }

  bbLimits[NUM_BARK_BANDS] = nSpectrumLength-1;

  return bbLimits;
};

Loudness.prototype.calcSpecificLoudness = function(bbLimits) {
  var NUM_BARK_BANDS = this.NUM_BARK_BANDS;
  var specific = this.specific;
  var normalisedSpectrum = this.normalisedSpectrum;

  for (var i = 0; i < NUM_BARK_BANDS; i++){
   var sum = 0;
   for (var j = bbLimits[i] ; j < bbLimits[i+1] ; j++) {

     sum += normalisedSpectrum[j];
   }

   specific[i] = Math.pow(sum, 0.23);
  }

  return specific;
};

Loudness.prototype.sumArray = function(array) {
  var sum = 0;

  for (var i = 0; i < array.length; i++) {
    sum += array[i];
  }

  return sum;
};

Loudness.prototype.process = function(buffer) {
    
  var bbLimits  = this.createBarkBandLimits();

  //process
  var specific  = this.calcSpecificLoudness(bbLimits);
  var total = this.sumArray(specific);

  return {
    specific: specific,
    total: total
  };
};

module.exports = Loudness;
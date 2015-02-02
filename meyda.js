!function(e){if("object"==typeof exports&&"undefined"!=typeof module)module.exports=e();else if("function"==typeof define&&define.amd)define([],e);else{var f;"undefined"!=typeof window?f=window:"undefined"!=typeof global?f=global:"undefined"!=typeof self&&(f=self),f.Meyda=e()}}(function(){var define,module,exports;return (function e(t,n,r){function s(o,u){if(!n[o]){if(!t[o]){var a=typeof require=="function"&&require;if(!u&&a)return a(o,!0);if(i)return i(o,!0);var f=new Error("Cannot find module '"+o+"'");throw f.code="MODULE_NOT_FOUND",f}var l=n[o]={exports:{}};t[o][0].call(l.exports,function(e){var n=t[o][1][e];return s(n?n:e)},l,l.exports,e,t,n,r)}return n[o].exports}var i=typeof require=="function"&&require;for(var o=0;o<r.length;o++)s(r[o]);return s})({1:[function(require,module,exports){
'use strict';

!function(exports, undefined) {

  var
    // If the typed array is unspecified, use this.
    DefaultArrayType = Float32Array,
    // Simple math functions we need.
    sqrt = Math.sqrt,
    sqr = function(number) {return Math.pow(number, 2)},
    // Internal convenience copies of the exported functions
    isComplexArray,
    ComplexArray

  exports.isComplexArray = isComplexArray = function(obj) {
    return obj !== undefined &&
      obj.hasOwnProperty !== undefined &&
      obj.hasOwnProperty('real') &&
      obj.hasOwnProperty('imag')
  }

  exports.ComplexArray = ComplexArray = function(other, opt_array_type){
    if (isComplexArray(other)) {
      // Copy constuctor.
      this.ArrayType = other.ArrayType
      this.real = new this.ArrayType(other.real)
      this.imag = new this.ArrayType(other.imag)
    } else {
      this.ArrayType = opt_array_type || DefaultArrayType
      // other can be either an array or a number.
      this.real = new this.ArrayType(other)
      this.imag = new this.ArrayType(this.real.length)
    }

    this.length = this.real.length
  }

  ComplexArray.prototype.toString = function() {
    var components = []

    this.forEach(function(c_value, i) {
      components.push(
        '(' +
        c_value.real.toFixed(2) + ',' +
        c_value.imag.toFixed(2) +
        ')'
      )
    })

    return '[' + components.join(',') + ']'
  }

  // In-place mapper.
  ComplexArray.prototype.map = function(mapper) {
    var
      i,
      n = this.length,
      // For GC efficiency, pass a single c_value object to the mapper.
      c_value = {}

    for (i = 0; i < n; i++) {
      c_value.real = this.real[i]
      c_value.imag = this.imag[i]
      mapper(c_value, i, n)
      this.real[i] = c_value.real
      this.imag[i] = c_value.imag
    }

    return this
  }

  ComplexArray.prototype.forEach = function(iterator) {
    var
      i,
      n = this.length,
      // For consistency with .map.
      c_value = {}

    for (i = 0; i < n; i++) {
      c_value.real = this.real[i]
      c_value.imag = this.imag[i]
      iterator(c_value, i, n)
    }
  }

  ComplexArray.prototype.conjugate = function() {
    return (new ComplexArray(this)).map(function(value) {
      value.imag *= -1
    })
  }

  // Helper so we can make ArrayType objects returned have similar interfaces
  //   to ComplexArrays.
  function iterable(obj) {
    if (!obj.forEach)
      obj.forEach = function(iterator) {
        var i, n = this.length

        for (i = 0; i < n; i++)
          iterator(this[i], i, n)
      }

    return obj
  }

  ComplexArray.prototype.magnitude = function() {
    var mags = new this.ArrayType(this.length)

    this.forEach(function(value, i) {
      mags[i] = sqrt(sqr(value.real) + sqr(value.imag))
    })

    // ArrayType will not necessarily be iterable: make it so.
    return iterable(mags)
  }
}(typeof exports === 'undefined' && (this.complex_array = {}) || exports)

},{}],2:[function(require,module,exports){
'use strict';

!function(exports, complex_array) {

  var
    ComplexArray = complex_array.ComplexArray,
    // Math constants and functions we need.
    PI = Math.PI,
    SQRT1_2 = Math.SQRT1_2,
    sqrt = Math.sqrt,
    cos = Math.cos,
    sin = Math.sin

  ComplexArray.prototype.FFT = function() {
    return FFT(this, false);
  }

  exports.FFT = function(input) {
    return ensureComplexArray(input).FFT()
  }

  ComplexArray.prototype.InvFFT = function() {
    return FFT(this, true)
  }

  exports.InvFFT = function(input) {
    return ensureComplexArray(input).InvFFT()
  }

  // Applies a frequency-space filter to input, and returns the real-space
  // filtered input.
  // filterer accepts freq, i, n and modifies freq.real and freq.imag.
  ComplexArray.prototype.frequencyMap = function(filterer) {
    return this.FFT().map(filterer).InvFFT()
  }

  exports.frequencyMap = function(input, filterer) {
    return ensureComplexArray(input).frequencyMap(filterer)
  }

  function ensureComplexArray(input) {
    return complex_array.isComplexArray(input) && input ||
        new ComplexArray(input)
  }

  function FFT(input, inverse) {
    var n = input.length

    if (n & (n - 1)) {
      return FFT_Recursive(input, inverse)
    } else {
      return FFT_2_Iterative(input, inverse)
    }
  }

  function FFT_Recursive(input, inverse) {
    var
      n = input.length,
      // Counters.
      i, j,
      output,
      // Complex multiplier and its delta.
      f_r, f_i, del_f_r, del_f_i,
      // Lowest divisor and remainder.
      p, m,
      normalisation,
      recursive_result,
      _swap, _real, _imag

    if (n === 1) {
      return input
    }

    output = new ComplexArray(n, input.ArrayType)

    // Use the lowest odd factor, so we are able to use FFT_2_Iterative in the
    // recursive transforms optimally.
    p = LowestOddFactor(n)
    m = n / p
    normalisation = 1 / sqrt(p)
    recursive_result = new ComplexArray(m, input.ArrayType)

    // Loops go like O(n Σ p_i), where p_i are the prime factors of n.
    // for a power of a prime, p, this reduces to O(n p log_p n)
    for(j = 0; j < p; j++) {
      for(i = 0; i < m; i++) {
        recursive_result.real[i] = input.real[i * p + j]
        recursive_result.imag[i] = input.imag[i * p + j]
      }
      // Don't go deeper unless necessary to save allocs.
      if (m > 1) {
        recursive_result = FFT(recursive_result, inverse)
      }

      del_f_r = cos(2*PI*j/n)
      del_f_i = (inverse ? -1 : 1) * sin(2*PI*j/n)
      f_r = 1
      f_i = 0

      for(i = 0; i < n; i++) {
        _real = recursive_result.real[i % m]
        _imag = recursive_result.imag[i % m]

        output.real[i] += f_r * _real - f_i * _imag
        output.imag[i] += f_r * _imag + f_i * _real

        _swap = f_r * del_f_r - f_i * del_f_i
        f_i = f_r * del_f_i + f_i * del_f_r
        f_r = _swap
      }
    }

    // Copy back to input to match FFT_2_Iterative in-placeness
    // TODO: faster way of making this in-place?
    for(i = 0; i < n; i++) {
      input.real[i] = normalisation * output.real[i]
      input.imag[i] = normalisation * output.imag[i]
    }

    return input
  }

  function FFT_2_Iterative(input, inverse) {
    var
      n = input.length,
      // Counters.
      i, j,
      output, output_r, output_i,
      // Complex multiplier and its delta.
      f_r, f_i, del_f_r, del_f_i, temp,
      // Temporary loop variables.
      l_index, r_index,
      left_r, left_i, right_r, right_i,
      // width of each sub-array for which we're iteratively calculating FFT.
      width

    output = BitReverseComplexArray(input)
    output_r = output.real
    output_i = output.imag
    // Loops go like O(n log n):
    //   width ~ log n; i,j ~ n
    width = 1
    while (width < n) {
      del_f_r = cos(PI/width)
      del_f_i = (inverse ? -1 : 1) * sin(PI/width)
      for (i = 0; i < n/(2*width); i++) {
        f_r = 1
        f_i = 0
        for (j = 0; j < width; j++) {
          l_index = 2*i*width + j
          r_index = l_index + width

          left_r = output_r[l_index]
          left_i = output_i[l_index]
          right_r = f_r * output_r[r_index] - f_i * output_i[r_index]
          right_i = f_i * output_r[r_index] + f_r * output_i[r_index]

          output_r[l_index] = SQRT1_2 * (left_r + right_r)
          output_i[l_index] = SQRT1_2 * (left_i + right_i)
          output_r[r_index] = SQRT1_2 * (left_r - right_r)
          output_i[r_index] = SQRT1_2 * (left_i - right_i)
          temp = f_r * del_f_r - f_i * del_f_i
          f_i = f_r * del_f_i + f_i * del_f_r
          f_r = temp
        }
      }
      width <<= 1
    }

    return output
  }

  function BitReverseIndex(index, n) {
    var bitreversed_index = 0

    while (n > 1) {
      bitreversed_index <<= 1
      bitreversed_index += index & 1
      index >>= 1
      n >>= 1
    }
    return bitreversed_index
  }

  function BitReverseComplexArray(array) {
    var n = array.length,
        flips = {},
        swap,
        i

    for(i = 0; i < n; i++) {
      var r_i = BitReverseIndex(i, n)

      if (flips.hasOwnProperty(i) || flips.hasOwnProperty(r_i)) continue

      swap = array.real[r_i]
      array.real[r_i] = array.real[i]
      array.real[i] = swap

      swap = array.imag[r_i]
      array.imag[r_i] = array.imag[i]
      array.imag[i] = swap

      flips[i] = flips[r_i] = true
    }

    return array
  }

  function LowestOddFactor(n) {
    var factor = 3,
        sqrt_n = sqrt(n)

    while(factor <= sqrt_n) {
      if (n % factor === 0) return factor
      factor = factor + 2
    }
    return n
  }

}(
  typeof exports === 'undefined' && (this.fft = {}) || exports,
  typeof require === 'undefined' && (this.complex_array) ||
    require('./complex_array')
)

},{"./complex_array":1}],3:[function(require,module,exports){
module.exports = {
  // buffer: require('./buffer'),
  // rms: require('./rms'),
  // energy: require('./energy'),
  // complexSpectrum: require('./complexSpectrum'),
  // spectralSlope: require('./spectralSlope'),
  // spectralCentroid: require('./spectralCentroid'),
  // spectralRolloff: require('./spectralRolloff'),
  // spectralFlatness: require('./spectralFlatness'),
  // spectralSpread: require('./spectralSpread'),
  // spectralSkewness: require('./spectralSkewness'),
  // spectralKurtosis: require('./spectralKurtosis'),
  // amplitudeSpectrum: require('./amplitudeSpectrum'),
  // zcr: require('./zcr'),
  // powerSpectrum: require('./powerSpectrum'),
  loudness: require('./loudness'),
  // perceptualSpread: require('./perceptualSpread'),
  // perceptualSharpness: require('./perceptualSharpness'),
  // mfcc: require('./mfcc')
};
},{"./loudness":4}],4:[function(require,module,exports){

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
},{}],5:[function(require,module,exports){
module.exportts = {

  "buffer": {
    "type": "array"
  },
  "rms": {
    "type": "number"
  },
  "energy": {
    "type": "number"
  },
  "zcr": {
    "type": "number"
  },
  "complexSpectrum": {
    "type": "multipleArrays",
    "arrayNames": {
      "1": "real",
      "2": "imag"
    }
  },
  "amplitudeSpectrum": {
    "type": "array"
  },
  "powerSpectrum": {
    "type": "array"
  },
  "spectralCentroid": {
    "type": "number"
  },
  "spectralFlatness": {
    "type": "number"
  },
  "spectralSlope": {
    "type": "number"
  },
  "spectralRolloff": {
    "type": "number"
  },
  "spectralSpread": {
    "type": "number"
  },
  "spectralSkewness": {
    "type": "number"
  },
  "spectralKurtosis": {
    "type": "number"
  },
  "loudness": {
    "type": "multipleArrays",
    "arrayNames": {
      "1": "total",
      "2": "specific"
    }
  },
  "perceptualSpread": {
    "type": "number"
  },
  "perceptualSharpness": {
    "type": "number"
  },
  "mfcc": {
    "type": "array"
  }
};
},{}],6:[function(require,module,exports){
// Meyda Javascript DSP library

// Dependencies
// ------------

var ComplexArray 	= require('../lib/jsfft/complex_array').ComplexArray,
		fft 					= require('../lib/jsfft/fft'); // modifies ComplexArray

var utils 				= require('./utils'),
		isPowerOfTwo 	= utils.isPowerOfTwo,
		µ 						= utils.µ;


// Constructor
// -----------

function Meyda(audioContext, src, bufSize, callback) {
	// reference to use inside new scopes
	var self = this;
	
	// make sure we can work
	if(!isPowerOfTwo(bufSize)) {
		throw new Error("Buffer size is not a power of two: Meyda will not run.");
	}
	
	if(!audioContext) {
		throw new Error("AudioContext wasn't specified: Meyda will not run.");
	}

	var bufferSize = bufSize || 256; //default buffer size
	var sampleRate = audioContext.sampleRate;
	var source = src; //initial source

	this.audioContext = audioContext;
	this.featureExtractors = {};

	//callback controllers
	this.EXTRACTION_STARTED = false;
	this._featuresToExtract = null;
	
	//WINDOWING
	//set default
	this.windowingFunction = "hanning";

	//initilize bark scale (to be used in most perceptual features).
	this.barkScale = this.computeBarkScale(bufferSize, sampleRate);

	//create windows
	this.hanning = this.computeHanning(bufferSize);
	this.hamming = this.computeHamming(bufferSize);
	// this.blackman = this.computeBlackman(bufferSize);

	this.featureInfo = require('./feature-info'); // to be included per module

	//create complexarray to hold the spectrum
	var computedSpectrumData = this.computeSpectrumData(bufferSize);

	//assign to meyda
	this.spectrumData = computedSpectrumData.data;
	this.complexSpectrum = computedSpectrumData.spectrum;
	this.ampSpectrum = computedSpectrumData.ampSpectrum;

	// console.log(computedSpectrumData);

		// console.log(ampSpectrum)
	this.initialiseExtractors();

	//create nodes
	window.spn = audioContext.createScriptProcessor(bufferSize,1,1);

	window.spn.onaudioprocess = function(e) {
		// "self" land over here because of the function closure

		//this is to obtain the current amplitude spectrum
		var signal = self.signal = e.inputBuffer.getChannelData(0);
		var data = self.spectrumData;
		var spec = self.complexSpectrum;
		var ampSpectrum = self.ampSpectrum;
		var windowedSignal = self.computeWindow(signal, self.windowingFunction);

		//map time domain
		data.map(function(value, i, n) {
			value.real = windowedSignal[i];
		});

		//calculate amplitude
		self.computeAmplitude(spec, ampSpectrum, bufferSize);

		//call callback if applicable
		if (typeof callback === "function" && EXTRACTION_STARTED) {
			callback(self.get(self._featuresToExtract));
		}

	};

	window.spn.connect(audioContext.destination);
	source.connect(window.spn, 0, 0);

	// constructors return "this" by default
}


// Compute methods
// ---------------

Meyda.prototype.computeAmplitude = function(complexSpectrum, ampSpectrum, bufferSize) {

	// works inside and outside Meyda
	complexSpectrum = complexSpectrum || this.complexSpectrum;
	ampSpectrum 		= ampSpectrum 		|| this.ampSpectrum;
	bufferSize 			= bufferSize 			|| this.bufferSize;

	for (var i = 0; i < bufferSize/2; i++) {
		ampSpectrum[i] = Math.sqrt(Math.pow(complexSpectrum.real[i],2) + Math.pow(complexSpectrum.imag[i],2));
	}
};

Meyda.prototype.computeHamming = function(bufferSize) {
	bufferSize = bufferSize || this.bufferSize;

	var hamming = new Float32Array(bufferSize);
	for (var i = 0; i < bufferSize; i++) {
		//According to http://uk.mathworks.com/help/signal/ref/hamming.html
		hamming[i] = 0.54 - 0.46*Math.cos(2*Math.PI*(i/bufferSize-1));
	}

	return hamming;
};

Meyda.prototype.computeHanning = function(bufferSize) {
	bufferSize = bufferSize || this.bufferSize;

	var hanning = new Float32Array(bufferSize);
	for (var i = 0; i < bufferSize; i++) {
		//According to the R documentation http://rgm.ogalab.net/RGM/R_rdfile?f=GENEAread/man/hanning.window.Rd&d=R_CC
		hanning[i] = 0.5 - 0.5*Math.cos(2*Math.PI*i/(bufferSize-1));
	}

	return hanning;
};

//UNFINISHED - blackman window implementation
/*
Meyda.prototype.computeBlackman = function(bufferSize) {
	bufferSize = bufferSize || this.bufferSize;
	
	var blackman = new Float32Array(bufferSize);
	//According to http://uk.mathworks.com/help/signal/ref/blackman.html
	//first half of the window
	for (var i = 0; i < (bufferSize % 2) ? (bufferSize+1)/2 : bufferSize/2; i++) {
		this.blackman[i] = 0.42 - 0.5*Math.cos(2*Math.PI*i/(bufferSize-1)) + 0.08*Math.cos(4*Math.PI*i/(bufferSize-1));
	}
	//second half of the window
	for (var i = bufferSize/2; i > 0; i--) {
		this.blackman[bufferSize - i] = this.blackman[i];
	}
};
*/

Meyda.prototype.computeWindow = function(sig, type) {

	var i, len = sig.length;
	var windowed = new Float32Array(len);

	for (i = 0; i < len; i++) {
		windowed[i] = sig[i] * this[type][i];
	}
	
	return windowed;
};

Meyda.prototype.computeBarkScale = function(bufferSize, sampleRate) {
	bufferSize = bufferSize || this.bufferSize;
	sampleRate = sampleRate || this.sampleRate;

  var barkScale = new Float32Array(bufferSize);

  for(var i = 0; i < bufferSize; i++){
    barkScale[i] = i * sampleRate / (bufferSize);
    barkScale[i] = 13 * Math.atan(barkScale[i]/1315.8) + 3.5 * Math.atan(Math.pow((barkScale[i]/7518), 2));
  }

  return barkScale;
};

Meyda.prototype.computeSpectrumData = function(bufferSize) {
	bufferSize = bufferSize || this.bufferSize;

	//create complexarray to hold the spectrum
	var data = new ComplexArray(bufferSize);
	var spectrum = data.FFT(); //transform
	var ampSpectrum = new Float32Array(bufferSize/2);

	return {
		data: data,
		spectrum: spectrum,
		ampSpectrum: ampSpectrum
	};
};


// Meyda methods
// -------------

// loads all the extractor objects
// initializes them and binds them
// to the featureExtractors
// and featureInfo lists
// @NOTE (c/sh)ould be handeled differently
Meyda.prototype.initialiseExtractors = function() {

	var extractors = require('./extractors');

	// Loudness
	var loudness = extractors.loudness({
		NUM_BARK_BANDS: 24,
		barkScale: this.barkScale,
		normalisedSpectrum: this.ampSpectrum,
		sampleRate: this.audioContext.sampleRate
	});

	this.featureExtractors.loudness = loudness;
	this.featureInfo.loudness = loudness.info;

	// Rest of the extractors
	// ...
};

//source setter method
Meyda.prototype.setSource = function(_src) {
	_src.connect(window.spn);
};

Meyda.prototype.start = function(features) {
	this._featuresToExtract = features;
	this.EXTRACTION_STARTED = true;
};

Meyda.prototype.stop = function() {
	this._featuresToExtract = null;
	this.EXTRACTION_STARTED = false;
};

// data pulling
Meyda.prototype.get = function(feature) {

	if(typeof feature === "object"){
		var results = {};
		for (var x = 0; x < feature.length; x++){
			try{
				results[feature[x]] = (this.featureExtractors[feature[x]].process(this.signal));
			} catch (e){
				console.error(e);
			}
		}
		return results;
	} else if (typeof feature === "string"){
		return this.featureExtractors[feature].process(this.signal);
	} else{
		throw new Error("Invalid Feature Format");
	}
};

module.exports = Meyda;

},{"../lib/jsfft/complex_array":1,"../lib/jsfft/fft":2,"./extractors":3,"./feature-info":5,"./utils":7}],7:[function(require,module,exports){
module.exports.µ = function(i, amplitudeSpect){
  var numerator = 0;
  var denominator = 0;
  
  for(var k = 0; k < amplitudeSpect.length; k++){
    numerator += Math.pow(k,i)*Math.abs(amplitudeSpect[k]);
    denominator += amplitudeSpect[k];
  }

  return numerator/denominator;
};

module.exports.isPowerOfTwo = function(num) {
  while (((num % 2) === 0) && num > 1) {
    num /= 2;
  }

  return (num == 1);
};
},{}]},{},[6])(6)
});
//# sourceMappingURL=data:application/json;base64,eyJ2ZXJzaW9uIjozLCJzb3VyY2VzIjpbIi4uLy4uLy4uLy4uL3Vzci9sb2NhbC9saWIvbm9kZV9tb2R1bGVzL2Jyb3dzZXJpZnkvbm9kZV9tb2R1bGVzL2Jyb3dzZXItcGFjay9fcHJlbHVkZS5qcyIsImxpYi9qc2ZmdC9jb21wbGV4X2FycmF5LmpzIiwibGliL2pzZmZ0L2ZmdC5qcyIsInNyYy9leHRyYWN0b3JzL2luZGV4LmpzIiwic3JjL2V4dHJhY3RvcnMvbG91ZG5lc3MuanMiLCJzcmMvZmVhdHVyZS1pbmZvLmpzIiwic3JjL21leWRhLmpzIiwic3JjL3V0aWxzLmpzIl0sIm5hbWVzIjpbXSwibWFwcGluZ3MiOiJBQUFBO0FDQUE7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBOztBQ3BIQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTs7QUNqT0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTs7QUNuQkE7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTs7QUNqR0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTs7QUNoRUE7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7O0FDeFFBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBIiwiZmlsZSI6ImdlbmVyYXRlZC5qcyIsInNvdXJjZVJvb3QiOiIiLCJzb3VyY2VzQ29udGVudCI6WyIoZnVuY3Rpb24gZSh0LG4scil7ZnVuY3Rpb24gcyhvLHUpe2lmKCFuW29dKXtpZighdFtvXSl7dmFyIGE9dHlwZW9mIHJlcXVpcmU9PVwiZnVuY3Rpb25cIiYmcmVxdWlyZTtpZighdSYmYSlyZXR1cm4gYShvLCEwKTtpZihpKXJldHVybiBpKG8sITApO3ZhciBmPW5ldyBFcnJvcihcIkNhbm5vdCBmaW5kIG1vZHVsZSAnXCIrbytcIidcIik7dGhyb3cgZi5jb2RlPVwiTU9EVUxFX05PVF9GT1VORFwiLGZ9dmFyIGw9bltvXT17ZXhwb3J0czp7fX07dFtvXVswXS5jYWxsKGwuZXhwb3J0cyxmdW5jdGlvbihlKXt2YXIgbj10W29dWzFdW2VdO3JldHVybiBzKG4/bjplKX0sbCxsLmV4cG9ydHMsZSx0LG4scil9cmV0dXJuIG5bb10uZXhwb3J0c312YXIgaT10eXBlb2YgcmVxdWlyZT09XCJmdW5jdGlvblwiJiZyZXF1aXJlO2Zvcih2YXIgbz0wO288ci5sZW5ndGg7bysrKXMocltvXSk7cmV0dXJuIHN9KSIsIid1c2Ugc3RyaWN0JztcblxuIWZ1bmN0aW9uKGV4cG9ydHMsIHVuZGVmaW5lZCkge1xuXG4gIHZhclxuICAgIC8vIElmIHRoZSB0eXBlZCBhcnJheSBpcyB1bnNwZWNpZmllZCwgdXNlIHRoaXMuXG4gICAgRGVmYXVsdEFycmF5VHlwZSA9IEZsb2F0MzJBcnJheSxcbiAgICAvLyBTaW1wbGUgbWF0aCBmdW5jdGlvbnMgd2UgbmVlZC5cbiAgICBzcXJ0ID0gTWF0aC5zcXJ0LFxuICAgIHNxciA9IGZ1bmN0aW9uKG51bWJlcikge3JldHVybiBNYXRoLnBvdyhudW1iZXIsIDIpfSxcbiAgICAvLyBJbnRlcm5hbCBjb252ZW5pZW5jZSBjb3BpZXMgb2YgdGhlIGV4cG9ydGVkIGZ1bmN0aW9uc1xuICAgIGlzQ29tcGxleEFycmF5LFxuICAgIENvbXBsZXhBcnJheVxuXG4gIGV4cG9ydHMuaXNDb21wbGV4QXJyYXkgPSBpc0NvbXBsZXhBcnJheSA9IGZ1bmN0aW9uKG9iaikge1xuICAgIHJldHVybiBvYmogIT09IHVuZGVmaW5lZCAmJlxuICAgICAgb2JqLmhhc093blByb3BlcnR5ICE9PSB1bmRlZmluZWQgJiZcbiAgICAgIG9iai5oYXNPd25Qcm9wZXJ0eSgncmVhbCcpICYmXG4gICAgICBvYmouaGFzT3duUHJvcGVydHkoJ2ltYWcnKVxuICB9XG5cbiAgZXhwb3J0cy5Db21wbGV4QXJyYXkgPSBDb21wbGV4QXJyYXkgPSBmdW5jdGlvbihvdGhlciwgb3B0X2FycmF5X3R5cGUpe1xuICAgIGlmIChpc0NvbXBsZXhBcnJheShvdGhlcikpIHtcbiAgICAgIC8vIENvcHkgY29uc3R1Y3Rvci5cbiAgICAgIHRoaXMuQXJyYXlUeXBlID0gb3RoZXIuQXJyYXlUeXBlXG4gICAgICB0aGlzLnJlYWwgPSBuZXcgdGhpcy5BcnJheVR5cGUob3RoZXIucmVhbClcbiAgICAgIHRoaXMuaW1hZyA9IG5ldyB0aGlzLkFycmF5VHlwZShvdGhlci5pbWFnKVxuICAgIH0gZWxzZSB7XG4gICAgICB0aGlzLkFycmF5VHlwZSA9IG9wdF9hcnJheV90eXBlIHx8IERlZmF1bHRBcnJheVR5cGVcbiAgICAgIC8vIG90aGVyIGNhbiBiZSBlaXRoZXIgYW4gYXJyYXkgb3IgYSBudW1iZXIuXG4gICAgICB0aGlzLnJlYWwgPSBuZXcgdGhpcy5BcnJheVR5cGUob3RoZXIpXG4gICAgICB0aGlzLmltYWcgPSBuZXcgdGhpcy5BcnJheVR5cGUodGhpcy5yZWFsLmxlbmd0aClcbiAgICB9XG5cbiAgICB0aGlzLmxlbmd0aCA9IHRoaXMucmVhbC5sZW5ndGhcbiAgfVxuXG4gIENvbXBsZXhBcnJheS5wcm90b3R5cGUudG9TdHJpbmcgPSBmdW5jdGlvbigpIHtcbiAgICB2YXIgY29tcG9uZW50cyA9IFtdXG5cbiAgICB0aGlzLmZvckVhY2goZnVuY3Rpb24oY192YWx1ZSwgaSkge1xuICAgICAgY29tcG9uZW50cy5wdXNoKFxuICAgICAgICAnKCcgK1xuICAgICAgICBjX3ZhbHVlLnJlYWwudG9GaXhlZCgyKSArICcsJyArXG4gICAgICAgIGNfdmFsdWUuaW1hZy50b0ZpeGVkKDIpICtcbiAgICAgICAgJyknXG4gICAgICApXG4gICAgfSlcblxuICAgIHJldHVybiAnWycgKyBjb21wb25lbnRzLmpvaW4oJywnKSArICddJ1xuICB9XG5cbiAgLy8gSW4tcGxhY2UgbWFwcGVyLlxuICBDb21wbGV4QXJyYXkucHJvdG90eXBlLm1hcCA9IGZ1bmN0aW9uKG1hcHBlcikge1xuICAgIHZhclxuICAgICAgaSxcbiAgICAgIG4gPSB0aGlzLmxlbmd0aCxcbiAgICAgIC8vIEZvciBHQyBlZmZpY2llbmN5LCBwYXNzIGEgc2luZ2xlIGNfdmFsdWUgb2JqZWN0IHRvIHRoZSBtYXBwZXIuXG4gICAgICBjX3ZhbHVlID0ge31cblxuICAgIGZvciAoaSA9IDA7IGkgPCBuOyBpKyspIHtcbiAgICAgIGNfdmFsdWUucmVhbCA9IHRoaXMucmVhbFtpXVxuICAgICAgY192YWx1ZS5pbWFnID0gdGhpcy5pbWFnW2ldXG4gICAgICBtYXBwZXIoY192YWx1ZSwgaSwgbilcbiAgICAgIHRoaXMucmVhbFtpXSA9IGNfdmFsdWUucmVhbFxuICAgICAgdGhpcy5pbWFnW2ldID0gY192YWx1ZS5pbWFnXG4gICAgfVxuXG4gICAgcmV0dXJuIHRoaXNcbiAgfVxuXG4gIENvbXBsZXhBcnJheS5wcm90b3R5cGUuZm9yRWFjaCA9IGZ1bmN0aW9uKGl0ZXJhdG9yKSB7XG4gICAgdmFyXG4gICAgICBpLFxuICAgICAgbiA9IHRoaXMubGVuZ3RoLFxuICAgICAgLy8gRm9yIGNvbnNpc3RlbmN5IHdpdGggLm1hcC5cbiAgICAgIGNfdmFsdWUgPSB7fVxuXG4gICAgZm9yIChpID0gMDsgaSA8IG47IGkrKykge1xuICAgICAgY192YWx1ZS5yZWFsID0gdGhpcy5yZWFsW2ldXG4gICAgICBjX3ZhbHVlLmltYWcgPSB0aGlzLmltYWdbaV1cbiAgICAgIGl0ZXJhdG9yKGNfdmFsdWUsIGksIG4pXG4gICAgfVxuICB9XG5cbiAgQ29tcGxleEFycmF5LnByb3RvdHlwZS5jb25qdWdhdGUgPSBmdW5jdGlvbigpIHtcbiAgICByZXR1cm4gKG5ldyBDb21wbGV4QXJyYXkodGhpcykpLm1hcChmdW5jdGlvbih2YWx1ZSkge1xuICAgICAgdmFsdWUuaW1hZyAqPSAtMVxuICAgIH0pXG4gIH1cblxuICAvLyBIZWxwZXIgc28gd2UgY2FuIG1ha2UgQXJyYXlUeXBlIG9iamVjdHMgcmV0dXJuZWQgaGF2ZSBzaW1pbGFyIGludGVyZmFjZXNcbiAgLy8gICB0byBDb21wbGV4QXJyYXlzLlxuICBmdW5jdGlvbiBpdGVyYWJsZShvYmopIHtcbiAgICBpZiAoIW9iai5mb3JFYWNoKVxuICAgICAgb2JqLmZvckVhY2ggPSBmdW5jdGlvbihpdGVyYXRvcikge1xuICAgICAgICB2YXIgaSwgbiA9IHRoaXMubGVuZ3RoXG5cbiAgICAgICAgZm9yIChpID0gMDsgaSA8IG47IGkrKylcbiAgICAgICAgICBpdGVyYXRvcih0aGlzW2ldLCBpLCBuKVxuICAgICAgfVxuXG4gICAgcmV0dXJuIG9ialxuICB9XG5cbiAgQ29tcGxleEFycmF5LnByb3RvdHlwZS5tYWduaXR1ZGUgPSBmdW5jdGlvbigpIHtcbiAgICB2YXIgbWFncyA9IG5ldyB0aGlzLkFycmF5VHlwZSh0aGlzLmxlbmd0aClcblxuICAgIHRoaXMuZm9yRWFjaChmdW5jdGlvbih2YWx1ZSwgaSkge1xuICAgICAgbWFnc1tpXSA9IHNxcnQoc3FyKHZhbHVlLnJlYWwpICsgc3FyKHZhbHVlLmltYWcpKVxuICAgIH0pXG5cbiAgICAvLyBBcnJheVR5cGUgd2lsbCBub3QgbmVjZXNzYXJpbHkgYmUgaXRlcmFibGU6IG1ha2UgaXQgc28uXG4gICAgcmV0dXJuIGl0ZXJhYmxlKG1hZ3MpXG4gIH1cbn0odHlwZW9mIGV4cG9ydHMgPT09ICd1bmRlZmluZWQnICYmICh0aGlzLmNvbXBsZXhfYXJyYXkgPSB7fSkgfHwgZXhwb3J0cylcbiIsIid1c2Ugc3RyaWN0JztcblxuIWZ1bmN0aW9uKGV4cG9ydHMsIGNvbXBsZXhfYXJyYXkpIHtcblxuICB2YXJcbiAgICBDb21wbGV4QXJyYXkgPSBjb21wbGV4X2FycmF5LkNvbXBsZXhBcnJheSxcbiAgICAvLyBNYXRoIGNvbnN0YW50cyBhbmQgZnVuY3Rpb25zIHdlIG5lZWQuXG4gICAgUEkgPSBNYXRoLlBJLFxuICAgIFNRUlQxXzIgPSBNYXRoLlNRUlQxXzIsXG4gICAgc3FydCA9IE1hdGguc3FydCxcbiAgICBjb3MgPSBNYXRoLmNvcyxcbiAgICBzaW4gPSBNYXRoLnNpblxuXG4gIENvbXBsZXhBcnJheS5wcm90b3R5cGUuRkZUID0gZnVuY3Rpb24oKSB7XG4gICAgcmV0dXJuIEZGVCh0aGlzLCBmYWxzZSk7XG4gIH1cblxuICBleHBvcnRzLkZGVCA9IGZ1bmN0aW9uKGlucHV0KSB7XG4gICAgcmV0dXJuIGVuc3VyZUNvbXBsZXhBcnJheShpbnB1dCkuRkZUKClcbiAgfVxuXG4gIENvbXBsZXhBcnJheS5wcm90b3R5cGUuSW52RkZUID0gZnVuY3Rpb24oKSB7XG4gICAgcmV0dXJuIEZGVCh0aGlzLCB0cnVlKVxuICB9XG5cbiAgZXhwb3J0cy5JbnZGRlQgPSBmdW5jdGlvbihpbnB1dCkge1xuICAgIHJldHVybiBlbnN1cmVDb21wbGV4QXJyYXkoaW5wdXQpLkludkZGVCgpXG4gIH1cblxuICAvLyBBcHBsaWVzIGEgZnJlcXVlbmN5LXNwYWNlIGZpbHRlciB0byBpbnB1dCwgYW5kIHJldHVybnMgdGhlIHJlYWwtc3BhY2VcbiAgLy8gZmlsdGVyZWQgaW5wdXQuXG4gIC8vIGZpbHRlcmVyIGFjY2VwdHMgZnJlcSwgaSwgbiBhbmQgbW9kaWZpZXMgZnJlcS5yZWFsIGFuZCBmcmVxLmltYWcuXG4gIENvbXBsZXhBcnJheS5wcm90b3R5cGUuZnJlcXVlbmN5TWFwID0gZnVuY3Rpb24oZmlsdGVyZXIpIHtcbiAgICByZXR1cm4gdGhpcy5GRlQoKS5tYXAoZmlsdGVyZXIpLkludkZGVCgpXG4gIH1cblxuICBleHBvcnRzLmZyZXF1ZW5jeU1hcCA9IGZ1bmN0aW9uKGlucHV0LCBmaWx0ZXJlcikge1xuICAgIHJldHVybiBlbnN1cmVDb21wbGV4QXJyYXkoaW5wdXQpLmZyZXF1ZW5jeU1hcChmaWx0ZXJlcilcbiAgfVxuXG4gIGZ1bmN0aW9uIGVuc3VyZUNvbXBsZXhBcnJheShpbnB1dCkge1xuICAgIHJldHVybiBjb21wbGV4X2FycmF5LmlzQ29tcGxleEFycmF5KGlucHV0KSAmJiBpbnB1dCB8fFxuICAgICAgICBuZXcgQ29tcGxleEFycmF5KGlucHV0KVxuICB9XG5cbiAgZnVuY3Rpb24gRkZUKGlucHV0LCBpbnZlcnNlKSB7XG4gICAgdmFyIG4gPSBpbnB1dC5sZW5ndGhcblxuICAgIGlmIChuICYgKG4gLSAxKSkge1xuICAgICAgcmV0dXJuIEZGVF9SZWN1cnNpdmUoaW5wdXQsIGludmVyc2UpXG4gICAgfSBlbHNlIHtcbiAgICAgIHJldHVybiBGRlRfMl9JdGVyYXRpdmUoaW5wdXQsIGludmVyc2UpXG4gICAgfVxuICB9XG5cbiAgZnVuY3Rpb24gRkZUX1JlY3Vyc2l2ZShpbnB1dCwgaW52ZXJzZSkge1xuICAgIHZhclxuICAgICAgbiA9IGlucHV0Lmxlbmd0aCxcbiAgICAgIC8vIENvdW50ZXJzLlxuICAgICAgaSwgaixcbiAgICAgIG91dHB1dCxcbiAgICAgIC8vIENvbXBsZXggbXVsdGlwbGllciBhbmQgaXRzIGRlbHRhLlxuICAgICAgZl9yLCBmX2ksIGRlbF9mX3IsIGRlbF9mX2ksXG4gICAgICAvLyBMb3dlc3QgZGl2aXNvciBhbmQgcmVtYWluZGVyLlxuICAgICAgcCwgbSxcbiAgICAgIG5vcm1hbGlzYXRpb24sXG4gICAgICByZWN1cnNpdmVfcmVzdWx0LFxuICAgICAgX3N3YXAsIF9yZWFsLCBfaW1hZ1xuXG4gICAgaWYgKG4gPT09IDEpIHtcbiAgICAgIHJldHVybiBpbnB1dFxuICAgIH1cblxuICAgIG91dHB1dCA9IG5ldyBDb21wbGV4QXJyYXkobiwgaW5wdXQuQXJyYXlUeXBlKVxuXG4gICAgLy8gVXNlIHRoZSBsb3dlc3Qgb2RkIGZhY3Rvciwgc28gd2UgYXJlIGFibGUgdG8gdXNlIEZGVF8yX0l0ZXJhdGl2ZSBpbiB0aGVcbiAgICAvLyByZWN1cnNpdmUgdHJhbnNmb3JtcyBvcHRpbWFsbHkuXG4gICAgcCA9IExvd2VzdE9kZEZhY3RvcihuKVxuICAgIG0gPSBuIC8gcFxuICAgIG5vcm1hbGlzYXRpb24gPSAxIC8gc3FydChwKVxuICAgIHJlY3Vyc2l2ZV9yZXN1bHQgPSBuZXcgQ29tcGxleEFycmF5KG0sIGlucHV0LkFycmF5VHlwZSlcblxuICAgIC8vIExvb3BzIGdvIGxpa2UgTyhuIM6jIHBfaSksIHdoZXJlIHBfaSBhcmUgdGhlIHByaW1lIGZhY3RvcnMgb2Ygbi5cbiAgICAvLyBmb3IgYSBwb3dlciBvZiBhIHByaW1lLCBwLCB0aGlzIHJlZHVjZXMgdG8gTyhuIHAgbG9nX3AgbilcbiAgICBmb3IoaiA9IDA7IGogPCBwOyBqKyspIHtcbiAgICAgIGZvcihpID0gMDsgaSA8IG07IGkrKykge1xuICAgICAgICByZWN1cnNpdmVfcmVzdWx0LnJlYWxbaV0gPSBpbnB1dC5yZWFsW2kgKiBwICsgal1cbiAgICAgICAgcmVjdXJzaXZlX3Jlc3VsdC5pbWFnW2ldID0gaW5wdXQuaW1hZ1tpICogcCArIGpdXG4gICAgICB9XG4gICAgICAvLyBEb24ndCBnbyBkZWVwZXIgdW5sZXNzIG5lY2Vzc2FyeSB0byBzYXZlIGFsbG9jcy5cbiAgICAgIGlmIChtID4gMSkge1xuICAgICAgICByZWN1cnNpdmVfcmVzdWx0ID0gRkZUKHJlY3Vyc2l2ZV9yZXN1bHQsIGludmVyc2UpXG4gICAgICB9XG5cbiAgICAgIGRlbF9mX3IgPSBjb3MoMipQSSpqL24pXG4gICAgICBkZWxfZl9pID0gKGludmVyc2UgPyAtMSA6IDEpICogc2luKDIqUEkqai9uKVxuICAgICAgZl9yID0gMVxuICAgICAgZl9pID0gMFxuXG4gICAgICBmb3IoaSA9IDA7IGkgPCBuOyBpKyspIHtcbiAgICAgICAgX3JlYWwgPSByZWN1cnNpdmVfcmVzdWx0LnJlYWxbaSAlIG1dXG4gICAgICAgIF9pbWFnID0gcmVjdXJzaXZlX3Jlc3VsdC5pbWFnW2kgJSBtXVxuXG4gICAgICAgIG91dHB1dC5yZWFsW2ldICs9IGZfciAqIF9yZWFsIC0gZl9pICogX2ltYWdcbiAgICAgICAgb3V0cHV0LmltYWdbaV0gKz0gZl9yICogX2ltYWcgKyBmX2kgKiBfcmVhbFxuXG4gICAgICAgIF9zd2FwID0gZl9yICogZGVsX2ZfciAtIGZfaSAqIGRlbF9mX2lcbiAgICAgICAgZl9pID0gZl9yICogZGVsX2ZfaSArIGZfaSAqIGRlbF9mX3JcbiAgICAgICAgZl9yID0gX3N3YXBcbiAgICAgIH1cbiAgICB9XG5cbiAgICAvLyBDb3B5IGJhY2sgdG8gaW5wdXQgdG8gbWF0Y2ggRkZUXzJfSXRlcmF0aXZlIGluLXBsYWNlbmVzc1xuICAgIC8vIFRPRE86IGZhc3RlciB3YXkgb2YgbWFraW5nIHRoaXMgaW4tcGxhY2U/XG4gICAgZm9yKGkgPSAwOyBpIDwgbjsgaSsrKSB7XG4gICAgICBpbnB1dC5yZWFsW2ldID0gbm9ybWFsaXNhdGlvbiAqIG91dHB1dC5yZWFsW2ldXG4gICAgICBpbnB1dC5pbWFnW2ldID0gbm9ybWFsaXNhdGlvbiAqIG91dHB1dC5pbWFnW2ldXG4gICAgfVxuXG4gICAgcmV0dXJuIGlucHV0XG4gIH1cblxuICBmdW5jdGlvbiBGRlRfMl9JdGVyYXRpdmUoaW5wdXQsIGludmVyc2UpIHtcbiAgICB2YXJcbiAgICAgIG4gPSBpbnB1dC5sZW5ndGgsXG4gICAgICAvLyBDb3VudGVycy5cbiAgICAgIGksIGosXG4gICAgICBvdXRwdXQsIG91dHB1dF9yLCBvdXRwdXRfaSxcbiAgICAgIC8vIENvbXBsZXggbXVsdGlwbGllciBhbmQgaXRzIGRlbHRhLlxuICAgICAgZl9yLCBmX2ksIGRlbF9mX3IsIGRlbF9mX2ksIHRlbXAsXG4gICAgICAvLyBUZW1wb3JhcnkgbG9vcCB2YXJpYWJsZXMuXG4gICAgICBsX2luZGV4LCByX2luZGV4LFxuICAgICAgbGVmdF9yLCBsZWZ0X2ksIHJpZ2h0X3IsIHJpZ2h0X2ksXG4gICAgICAvLyB3aWR0aCBvZiBlYWNoIHN1Yi1hcnJheSBmb3Igd2hpY2ggd2UncmUgaXRlcmF0aXZlbHkgY2FsY3VsYXRpbmcgRkZULlxuICAgICAgd2lkdGhcblxuICAgIG91dHB1dCA9IEJpdFJldmVyc2VDb21wbGV4QXJyYXkoaW5wdXQpXG4gICAgb3V0cHV0X3IgPSBvdXRwdXQucmVhbFxuICAgIG91dHB1dF9pID0gb3V0cHV0LmltYWdcbiAgICAvLyBMb29wcyBnbyBsaWtlIE8obiBsb2cgbik6XG4gICAgLy8gICB3aWR0aCB+IGxvZyBuOyBpLGogfiBuXG4gICAgd2lkdGggPSAxXG4gICAgd2hpbGUgKHdpZHRoIDwgbikge1xuICAgICAgZGVsX2ZfciA9IGNvcyhQSS93aWR0aClcbiAgICAgIGRlbF9mX2kgPSAoaW52ZXJzZSA/IC0xIDogMSkgKiBzaW4oUEkvd2lkdGgpXG4gICAgICBmb3IgKGkgPSAwOyBpIDwgbi8oMip3aWR0aCk7IGkrKykge1xuICAgICAgICBmX3IgPSAxXG4gICAgICAgIGZfaSA9IDBcbiAgICAgICAgZm9yIChqID0gMDsgaiA8IHdpZHRoOyBqKyspIHtcbiAgICAgICAgICBsX2luZGV4ID0gMippKndpZHRoICsgalxuICAgICAgICAgIHJfaW5kZXggPSBsX2luZGV4ICsgd2lkdGhcblxuICAgICAgICAgIGxlZnRfciA9IG91dHB1dF9yW2xfaW5kZXhdXG4gICAgICAgICAgbGVmdF9pID0gb3V0cHV0X2lbbF9pbmRleF1cbiAgICAgICAgICByaWdodF9yID0gZl9yICogb3V0cHV0X3Jbcl9pbmRleF0gLSBmX2kgKiBvdXRwdXRfaVtyX2luZGV4XVxuICAgICAgICAgIHJpZ2h0X2kgPSBmX2kgKiBvdXRwdXRfcltyX2luZGV4XSArIGZfciAqIG91dHB1dF9pW3JfaW5kZXhdXG5cbiAgICAgICAgICBvdXRwdXRfcltsX2luZGV4XSA9IFNRUlQxXzIgKiAobGVmdF9yICsgcmlnaHRfcilcbiAgICAgICAgICBvdXRwdXRfaVtsX2luZGV4XSA9IFNRUlQxXzIgKiAobGVmdF9pICsgcmlnaHRfaSlcbiAgICAgICAgICBvdXRwdXRfcltyX2luZGV4XSA9IFNRUlQxXzIgKiAobGVmdF9yIC0gcmlnaHRfcilcbiAgICAgICAgICBvdXRwdXRfaVtyX2luZGV4XSA9IFNRUlQxXzIgKiAobGVmdF9pIC0gcmlnaHRfaSlcbiAgICAgICAgICB0ZW1wID0gZl9yICogZGVsX2ZfciAtIGZfaSAqIGRlbF9mX2lcbiAgICAgICAgICBmX2kgPSBmX3IgKiBkZWxfZl9pICsgZl9pICogZGVsX2ZfclxuICAgICAgICAgIGZfciA9IHRlbXBcbiAgICAgICAgfVxuICAgICAgfVxuICAgICAgd2lkdGggPDw9IDFcbiAgICB9XG5cbiAgICByZXR1cm4gb3V0cHV0XG4gIH1cblxuICBmdW5jdGlvbiBCaXRSZXZlcnNlSW5kZXgoaW5kZXgsIG4pIHtcbiAgICB2YXIgYml0cmV2ZXJzZWRfaW5kZXggPSAwXG5cbiAgICB3aGlsZSAobiA+IDEpIHtcbiAgICAgIGJpdHJldmVyc2VkX2luZGV4IDw8PSAxXG4gICAgICBiaXRyZXZlcnNlZF9pbmRleCArPSBpbmRleCAmIDFcbiAgICAgIGluZGV4ID4+PSAxXG4gICAgICBuID4+PSAxXG4gICAgfVxuICAgIHJldHVybiBiaXRyZXZlcnNlZF9pbmRleFxuICB9XG5cbiAgZnVuY3Rpb24gQml0UmV2ZXJzZUNvbXBsZXhBcnJheShhcnJheSkge1xuICAgIHZhciBuID0gYXJyYXkubGVuZ3RoLFxuICAgICAgICBmbGlwcyA9IHt9LFxuICAgICAgICBzd2FwLFxuICAgICAgICBpXG5cbiAgICBmb3IoaSA9IDA7IGkgPCBuOyBpKyspIHtcbiAgICAgIHZhciByX2kgPSBCaXRSZXZlcnNlSW5kZXgoaSwgbilcblxuICAgICAgaWYgKGZsaXBzLmhhc093blByb3BlcnR5KGkpIHx8IGZsaXBzLmhhc093blByb3BlcnR5KHJfaSkpIGNvbnRpbnVlXG5cbiAgICAgIHN3YXAgPSBhcnJheS5yZWFsW3JfaV1cbiAgICAgIGFycmF5LnJlYWxbcl9pXSA9IGFycmF5LnJlYWxbaV1cbiAgICAgIGFycmF5LnJlYWxbaV0gPSBzd2FwXG5cbiAgICAgIHN3YXAgPSBhcnJheS5pbWFnW3JfaV1cbiAgICAgIGFycmF5LmltYWdbcl9pXSA9IGFycmF5LmltYWdbaV1cbiAgICAgIGFycmF5LmltYWdbaV0gPSBzd2FwXG5cbiAgICAgIGZsaXBzW2ldID0gZmxpcHNbcl9pXSA9IHRydWVcbiAgICB9XG5cbiAgICByZXR1cm4gYXJyYXlcbiAgfVxuXG4gIGZ1bmN0aW9uIExvd2VzdE9kZEZhY3RvcihuKSB7XG4gICAgdmFyIGZhY3RvciA9IDMsXG4gICAgICAgIHNxcnRfbiA9IHNxcnQobilcblxuICAgIHdoaWxlKGZhY3RvciA8PSBzcXJ0X24pIHtcbiAgICAgIGlmIChuICUgZmFjdG9yID09PSAwKSByZXR1cm4gZmFjdG9yXG4gICAgICBmYWN0b3IgPSBmYWN0b3IgKyAyXG4gICAgfVxuICAgIHJldHVybiBuXG4gIH1cblxufShcbiAgdHlwZW9mIGV4cG9ydHMgPT09ICd1bmRlZmluZWQnICYmICh0aGlzLmZmdCA9IHt9KSB8fCBleHBvcnRzLFxuICB0eXBlb2YgcmVxdWlyZSA9PT0gJ3VuZGVmaW5lZCcgJiYgKHRoaXMuY29tcGxleF9hcnJheSkgfHxcbiAgICByZXF1aXJlKCcuL2NvbXBsZXhfYXJyYXknKVxuKVxuIiwibW9kdWxlLmV4cG9ydHMgPSB7XG4gIC8vIGJ1ZmZlcjogcmVxdWlyZSgnLi9idWZmZXInKSxcbiAgLy8gcm1zOiByZXF1aXJlKCcuL3JtcycpLFxuICAvLyBlbmVyZ3k6IHJlcXVpcmUoJy4vZW5lcmd5JyksXG4gIC8vIGNvbXBsZXhTcGVjdHJ1bTogcmVxdWlyZSgnLi9jb21wbGV4U3BlY3RydW0nKSxcbiAgLy8gc3BlY3RyYWxTbG9wZTogcmVxdWlyZSgnLi9zcGVjdHJhbFNsb3BlJyksXG4gIC8vIHNwZWN0cmFsQ2VudHJvaWQ6IHJlcXVpcmUoJy4vc3BlY3RyYWxDZW50cm9pZCcpLFxuICAvLyBzcGVjdHJhbFJvbGxvZmY6IHJlcXVpcmUoJy4vc3BlY3RyYWxSb2xsb2ZmJyksXG4gIC8vIHNwZWN0cmFsRmxhdG5lc3M6IHJlcXVpcmUoJy4vc3BlY3RyYWxGbGF0bmVzcycpLFxuICAvLyBzcGVjdHJhbFNwcmVhZDogcmVxdWlyZSgnLi9zcGVjdHJhbFNwcmVhZCcpLFxuICAvLyBzcGVjdHJhbFNrZXduZXNzOiByZXF1aXJlKCcuL3NwZWN0cmFsU2tld25lc3MnKSxcbiAgLy8gc3BlY3RyYWxLdXJ0b3NpczogcmVxdWlyZSgnLi9zcGVjdHJhbEt1cnRvc2lzJyksXG4gIC8vIGFtcGxpdHVkZVNwZWN0cnVtOiByZXF1aXJlKCcuL2FtcGxpdHVkZVNwZWN0cnVtJyksXG4gIC8vIHpjcjogcmVxdWlyZSgnLi96Y3InKSxcbiAgLy8gcG93ZXJTcGVjdHJ1bTogcmVxdWlyZSgnLi9wb3dlclNwZWN0cnVtJyksXG4gIGxvdWRuZXNzOiByZXF1aXJlKCcuL2xvdWRuZXNzJyksXG4gIC8vIHBlcmNlcHR1YWxTcHJlYWQ6IHJlcXVpcmUoJy4vcGVyY2VwdHVhbFNwcmVhZCcpLFxuICAvLyBwZXJjZXB0dWFsU2hhcnBuZXNzOiByZXF1aXJlKCcuL3BlcmNlcHR1YWxTaGFycG5lc3MnKSxcbiAgLy8gbWZjYzogcmVxdWlyZSgnLi9tZmNjJylcbn07IiwiXG5mdW5jdGlvbiBMb3VkbmVzcyhvcHRzKSB7XG4gIGlmICghKHRoaXMgaW5zdGFuY2VvZiBMb3VkbmVzcykpIHJldHVybiBuZXcgTG91ZG5lc3Mob3B0cyk7XG5cbiAgLy8gZmVhdHVyZUluZm8gZm9yIGxvdWRuZXNzXG4gIHRoaXMuaW5mbyA9IHtcbiAgICBcInR5cGVcIjogXCJtdWx0aXBsZUFycmF5c1wiLFxuICAgIFwiYXJyYXlOYW1lc1wiOiB7XG4gICAgICBcIjFcIjogXCJ0b3RhbFwiLFxuICAgICAgXCIyXCI6IFwic3BlY2lmaWNcIlxuICAgIH1cbiAgfTtcblxuICB0aGlzLk5VTV9CQVJLX0JBTkRTID0gb3B0cy5OVU1fQkFSS19CQU5EUyB8fCAyNDtcbiAgdGhpcy5ub3JtYWxpc2VkU3BlY3RydW0gPSBvcHRzLm5vcm1hbGlzZWRTcGVjdHJ1bTtcbiAgdGhpcy5zYW1wbGVSYXRlID0gb3B0cy5zYW1wbGVSYXRlO1xuXG4gIHRoaXMuc3BlY2lmaWMgPSBuZXcgRmxvYXQzMkFycmF5KHRoaXMuTlVNX0JBUktfQkFORFMpO1xuICB0aGlzLnByb2Nlc3MgPSB0aGlzLnByb2Nlc3MuYmluZCh0aGlzKTsgLy8gb3B0aW1pemUgbGF0ZXJcbiAgdGhpcy5iYXJrU2NhbGUgPSBvcHRzLmJhcmtTY2FsZTtcbiAgdGhpcy5iYkxpbWl0cyAgPSB0aGlzLmNvbXB1dGVCYXJrQmFuZExpbWl0cyh0aGlzLmJhcmtTY2FsZSwgdGhpcy5ub3JtYWxpc2VkU3BlY3RydW0ubGVuZ3RoLCB0aGlzLk5VTV9CQVJLX0JBTkRTKTtcbn1cblxuTG91ZG5lc3MucHJvdG90eXBlLmNvbXB1dGVCYXJrQmFuZExpbWl0cyA9IGZ1bmN0aW9uKGJhcmtTY2FsZSwgblNwZWN0cnVtTGVuZ3RoLCBudW1fYmFya19iYW5kcykge1xuICBiYXJrU2NhbGUgICAgICAgPSBiYXJrU2NhbGUgICAgICAgfHwgdGhpcy5iYXJrU2NhbGU7XG4gIG5TcGVjdHJ1bUxlbmd0aCA9IG5TcGVjdHJ1bUxlbmd0aCB8fCB0aGlzLm5vcm1hbGlzZWRTcGVjdHJ1bS5sZW5ndGg7XG4gIG51bV9iYXJrX2JhbmRzICA9IG51bV9iYXJrX2JhbmRzICB8fCB0aGlzLk5VTV9CQVJLX0JBTkRTO1xuICBcbiAgdmFyIGN1cnJlbnRCYW5kRW5kID0gYmFya1NjYWxlW25TcGVjdHJ1bUxlbmd0aC0xXS9udW1fYmFya19iYW5kcztcbiAgdmFyIGN1cnJlbnRCYW5kID0gMTtcblxuICB2YXIgYmJMaW1pdHMgPSBuZXcgSW50MzJBcnJheShudW1fYmFya19iYW5kcysxKTtcbiAgYmJMaW1pdHNbMF0gPSAwO1xuXG4gIGZvcih2YXIgaSA9IDA7IGk8blNwZWN0cnVtTGVuZ3RoOyBpKyspe1xuICAgd2hpbGUoYmFya1NjYWxlW2ldID4gY3VycmVudEJhbmRFbmQpIHtcbiAgICAgYmJMaW1pdHNbY3VycmVudEJhbmQrK10gPSBpO1xuICAgICBjdXJyZW50QmFuZEVuZCA9IGN1cnJlbnRCYW5kKmJhcmtTY2FsZVtuU3BlY3RydW1MZW5ndGgtMV0vbnVtX2JhcmtfYmFuZHM7XG4gICB9XG4gIH1cblxuICBiYkxpbWl0c1tudW1fYmFya19iYW5kc10gPSBuU3BlY3RydW1MZW5ndGgtMTtcblxuICByZXR1cm4gYmJMaW1pdHM7XG59O1xuXG5Mb3VkbmVzcy5wcm90b3R5cGUuY29tcHV0ZVNwZWNpZmljTG91ZG5lc3MgPSBmdW5jdGlvbihiYkxpbWl0cywgc3BlY2lmaWMsIG5TcGVjdHJ1bSwgbnVtX2JhcmtfYmFuZHMpIHtcbiAgYmJMaW1pdHMgICAgICAgID0gYmJMaW1pdHMgICAgICAgfHzCoHRoaXMuYmJMaW1pdHM7XG4gIHNwZWNpZmljICAgICAgICA9IHNwZWNpZmljICAgICAgIHx8wqB0aGlzLnNwZWNpZmljO1xuICBuU3BlY3RydW0gICAgICAgPSBuU3BlY3RydW0gICAgICB8fCB0aGlzLm5vcm1hbGlzZWRTcGVjdHJ1bTtcbiAgbnVtX2JhcmtfYmFuZHMgID0gbnVtX2JhcmtfYmFuZHMgfHwgdGhpcy5OVU1fQkFSS19CQU5EUztcblxuICAvLyBjb25zb2xlLmxvZyhiYkxpbWl0cywgc3BlY2lmaWMsIG5TcGVjdHJ1bSwgbnVtX2JhcmtfYmFuZHMpO1xuXG4gIGZvciAodmFyIGkgPSAwOyBpIDwgbnVtX2JhcmtfYmFuZHM7IGkrKyl7XG4gICAgdmFyIHN1bSA9IDA7XG5cbiAgICBmb3IgKHZhciBqID0gYmJMaW1pdHNbaV0gOyBqIDwgYmJMaW1pdHNbaSsxXSA7IGorKykge1xuICAgICAgc3VtICs9IG5TcGVjdHJ1bVtqXTtcbiAgICB9XG5cbiAgICBzcGVjaWZpY1tpXSA9IE1hdGgucG93KHN1bSwgMC4yMyk7XG4gIH1cblxuICByZXR1cm4gc3BlY2lmaWM7XG59O1xuXG4vLyBiZWxvbmdzIGluIHV0aWxzXG5Mb3VkbmVzcy5wcm90b3R5cGUuc3VtQXJyYXkgPSBmdW5jdGlvbihhcnJheSkge1xuICB2YXIgc3VtID0gMDtcblxuICBmb3IgKHZhciBpID0gMDsgaSA8IGFycmF5Lmxlbmd0aDsgaSsrKSB7XG4gICAgc3VtICs9IGFycmF5W2ldO1xuICB9XG5cbiAgcmV0dXJuIHN1bTtcbn07XG5cbkxvdWRuZXNzLnByb3RvdHlwZS5wcm9jZXNzID0gZnVuY3Rpb24oKSB7XG4gIHZhciBiYXJrU2NhbGUgPSB0aGlzLmJhcmtTY2FsZTtcbiAgdmFyIG5TcGVjdHJ1bSA9IHRoaXMubm9ybWFsaXNlZFNwZWN0cnVtO1xuICB2YXIgblNwZWN0cnVtTGVuZ3RoID0gblNwZWN0cnVtLmxlbmd0aDtcbiAgdmFyIE5VTV9CQVJLX0JBTkRTID0gdGhpcy5OVU1fQkFSS19CQU5EUztcbiAgdmFyIHNwZWNpZmljID0gdGhpcy5zcGVjaWZpYztcbiAgdmFyIGJiTGltaXRzICA9IHRoaXMuYmJMaW1pdHM7XG5cbiAgLy9wcm9jZXNzXG4gIHZhciBzcGVjICA9IHRoaXMuY29tcHV0ZVNwZWNpZmljTG91ZG5lc3MoYmJMaW1pdHMsIHNwZWNpZmljLCBuU3BlY3RydW0sIE5VTV9CQVJLX0JBTkRTKTtcbiAgLy8gY29uc29sZS5sb2coc3BlYywgdGhpcy5zcGVjaWZpYyk7XG4gIHZhciB0b3RhbCA9IHRoaXMuc3VtQXJyYXkoc3BlYyk7XG5cbiAgcmV0dXJuIHtcbiAgICBzcGVjaWZpYzogc3BlYyxcbiAgICB0b3RhbDogdG90YWxcbiAgfTtcbn07XG5cbm1vZHVsZS5leHBvcnRzID0gTG91ZG5lc3M7IiwibW9kdWxlLmV4cG9ydHRzID0ge1xuXG4gIFwiYnVmZmVyXCI6IHtcbiAgICBcInR5cGVcIjogXCJhcnJheVwiXG4gIH0sXG4gIFwicm1zXCI6IHtcbiAgICBcInR5cGVcIjogXCJudW1iZXJcIlxuICB9LFxuICBcImVuZXJneVwiOiB7XG4gICAgXCJ0eXBlXCI6IFwibnVtYmVyXCJcbiAgfSxcbiAgXCJ6Y3JcIjoge1xuICAgIFwidHlwZVwiOiBcIm51bWJlclwiXG4gIH0sXG4gIFwiY29tcGxleFNwZWN0cnVtXCI6IHtcbiAgICBcInR5cGVcIjogXCJtdWx0aXBsZUFycmF5c1wiLFxuICAgIFwiYXJyYXlOYW1lc1wiOiB7XG4gICAgICBcIjFcIjogXCJyZWFsXCIsXG4gICAgICBcIjJcIjogXCJpbWFnXCJcbiAgICB9XG4gIH0sXG4gIFwiYW1wbGl0dWRlU3BlY3RydW1cIjoge1xuICAgIFwidHlwZVwiOiBcImFycmF5XCJcbiAgfSxcbiAgXCJwb3dlclNwZWN0cnVtXCI6IHtcbiAgICBcInR5cGVcIjogXCJhcnJheVwiXG4gIH0sXG4gIFwic3BlY3RyYWxDZW50cm9pZFwiOiB7XG4gICAgXCJ0eXBlXCI6IFwibnVtYmVyXCJcbiAgfSxcbiAgXCJzcGVjdHJhbEZsYXRuZXNzXCI6IHtcbiAgICBcInR5cGVcIjogXCJudW1iZXJcIlxuICB9LFxuICBcInNwZWN0cmFsU2xvcGVcIjoge1xuICAgIFwidHlwZVwiOiBcIm51bWJlclwiXG4gIH0sXG4gIFwic3BlY3RyYWxSb2xsb2ZmXCI6IHtcbiAgICBcInR5cGVcIjogXCJudW1iZXJcIlxuICB9LFxuICBcInNwZWN0cmFsU3ByZWFkXCI6IHtcbiAgICBcInR5cGVcIjogXCJudW1iZXJcIlxuICB9LFxuICBcInNwZWN0cmFsU2tld25lc3NcIjoge1xuICAgIFwidHlwZVwiOiBcIm51bWJlclwiXG4gIH0sXG4gIFwic3BlY3RyYWxLdXJ0b3Npc1wiOiB7XG4gICAgXCJ0eXBlXCI6IFwibnVtYmVyXCJcbiAgfSxcbiAgXCJsb3VkbmVzc1wiOiB7XG4gICAgXCJ0eXBlXCI6IFwibXVsdGlwbGVBcnJheXNcIixcbiAgICBcImFycmF5TmFtZXNcIjoge1xuICAgICAgXCIxXCI6IFwidG90YWxcIixcbiAgICAgIFwiMlwiOiBcInNwZWNpZmljXCJcbiAgICB9XG4gIH0sXG4gIFwicGVyY2VwdHVhbFNwcmVhZFwiOiB7XG4gICAgXCJ0eXBlXCI6IFwibnVtYmVyXCJcbiAgfSxcbiAgXCJwZXJjZXB0dWFsU2hhcnBuZXNzXCI6IHtcbiAgICBcInR5cGVcIjogXCJudW1iZXJcIlxuICB9LFxuICBcIm1mY2NcIjoge1xuICAgIFwidHlwZVwiOiBcImFycmF5XCJcbiAgfVxufTsiLCIvLyBNZXlkYSBKYXZhc2NyaXB0IERTUCBsaWJyYXJ5XG5cbi8vIERlcGVuZGVuY2llc1xuLy8gLS0tLS0tLS0tLS0tXG5cbnZhciBDb21wbGV4QXJyYXkgXHQ9IHJlcXVpcmUoJy4uL2xpYi9qc2ZmdC9jb21wbGV4X2FycmF5JykuQ29tcGxleEFycmF5LFxuXHRcdGZmdCBcdFx0XHRcdFx0PSByZXF1aXJlKCcuLi9saWIvanNmZnQvZmZ0Jyk7IC8vIG1vZGlmaWVzIENvbXBsZXhBcnJheVxuXG52YXIgdXRpbHMgXHRcdFx0XHQ9IHJlcXVpcmUoJy4vdXRpbHMnKSxcblx0XHRpc1Bvd2VyT2ZUd28gXHQ9IHV0aWxzLmlzUG93ZXJPZlR3byxcblx0XHTCtSBcdFx0XHRcdFx0XHQ9IHV0aWxzLsK1O1xuXG5cbi8vIENvbnN0cnVjdG9yXG4vLyAtLS0tLS0tLS0tLVxuXG5mdW5jdGlvbiBNZXlkYShhdWRpb0NvbnRleHQsIHNyYywgYnVmU2l6ZSwgY2FsbGJhY2spIHtcblx0Ly8gcmVmZXJlbmNlIHRvIHVzZSBpbnNpZGUgbmV3IHNjb3Blc1xuXHR2YXIgc2VsZiA9IHRoaXM7XG5cdFxuXHQvLyBtYWtlIHN1cmUgd2UgY2FuIHdvcmtcblx0aWYoIWlzUG93ZXJPZlR3byhidWZTaXplKSkge1xuXHRcdHRocm93IG5ldyBFcnJvcihcIkJ1ZmZlciBzaXplIGlzIG5vdCBhIHBvd2VyIG9mIHR3bzogTWV5ZGEgd2lsbCBub3QgcnVuLlwiKTtcblx0fVxuXHRcblx0aWYoIWF1ZGlvQ29udGV4dCkge1xuXHRcdHRocm93IG5ldyBFcnJvcihcIkF1ZGlvQ29udGV4dCB3YXNuJ3Qgc3BlY2lmaWVkOiBNZXlkYSB3aWxsIG5vdCBydW4uXCIpO1xuXHR9XG5cblx0dmFyIGJ1ZmZlclNpemUgPSBidWZTaXplIHx8IDI1NjsgLy9kZWZhdWx0IGJ1ZmZlciBzaXplXG5cdHZhciBzYW1wbGVSYXRlID0gYXVkaW9Db250ZXh0LnNhbXBsZVJhdGU7XG5cdHZhciBzb3VyY2UgPSBzcmM7IC8vaW5pdGlhbCBzb3VyY2VcblxuXHR0aGlzLmF1ZGlvQ29udGV4dCA9IGF1ZGlvQ29udGV4dDtcblx0dGhpcy5mZWF0dXJlRXh0cmFjdG9ycyA9IHt9O1xuXG5cdC8vY2FsbGJhY2sgY29udHJvbGxlcnNcblx0dGhpcy5FWFRSQUNUSU9OX1NUQVJURUQgPSBmYWxzZTtcblx0dGhpcy5fZmVhdHVyZXNUb0V4dHJhY3QgPSBudWxsO1xuXHRcblx0Ly9XSU5ET1dJTkdcblx0Ly9zZXQgZGVmYXVsdFxuXHR0aGlzLndpbmRvd2luZ0Z1bmN0aW9uID0gXCJoYW5uaW5nXCI7XG5cblx0Ly9pbml0aWxpemUgYmFyayBzY2FsZSAodG8gYmUgdXNlZCBpbiBtb3N0IHBlcmNlcHR1YWwgZmVhdHVyZXMpLlxuXHR0aGlzLmJhcmtTY2FsZSA9IHRoaXMuY29tcHV0ZUJhcmtTY2FsZShidWZmZXJTaXplLCBzYW1wbGVSYXRlKTtcblxuXHQvL2NyZWF0ZSB3aW5kb3dzXG5cdHRoaXMuaGFubmluZyA9IHRoaXMuY29tcHV0ZUhhbm5pbmcoYnVmZmVyU2l6ZSk7XG5cdHRoaXMuaGFtbWluZyA9IHRoaXMuY29tcHV0ZUhhbW1pbmcoYnVmZmVyU2l6ZSk7XG5cdC8vIHRoaXMuYmxhY2ttYW4gPSB0aGlzLmNvbXB1dGVCbGFja21hbihidWZmZXJTaXplKTtcblxuXHR0aGlzLmZlYXR1cmVJbmZvID0gcmVxdWlyZSgnLi9mZWF0dXJlLWluZm8nKTsgLy8gdG8gYmUgaW5jbHVkZWQgcGVyIG1vZHVsZVxuXG5cdC8vY3JlYXRlIGNvbXBsZXhhcnJheSB0byBob2xkIHRoZSBzcGVjdHJ1bVxuXHR2YXIgY29tcHV0ZWRTcGVjdHJ1bURhdGEgPSB0aGlzLmNvbXB1dGVTcGVjdHJ1bURhdGEoYnVmZmVyU2l6ZSk7XG5cblx0Ly9hc3NpZ24gdG8gbWV5ZGFcblx0dGhpcy5zcGVjdHJ1bURhdGEgPSBjb21wdXRlZFNwZWN0cnVtRGF0YS5kYXRhO1xuXHR0aGlzLmNvbXBsZXhTcGVjdHJ1bSA9IGNvbXB1dGVkU3BlY3RydW1EYXRhLnNwZWN0cnVtO1xuXHR0aGlzLmFtcFNwZWN0cnVtID0gY29tcHV0ZWRTcGVjdHJ1bURhdGEuYW1wU3BlY3RydW07XG5cblx0Ly8gY29uc29sZS5sb2coY29tcHV0ZWRTcGVjdHJ1bURhdGEpO1xuXG5cdFx0Ly8gY29uc29sZS5sb2coYW1wU3BlY3RydW0pXG5cdHRoaXMuaW5pdGlhbGlzZUV4dHJhY3RvcnMoKTtcblxuXHQvL2NyZWF0ZSBub2Rlc1xuXHR3aW5kb3cuc3BuID0gYXVkaW9Db250ZXh0LmNyZWF0ZVNjcmlwdFByb2Nlc3NvcihidWZmZXJTaXplLDEsMSk7XG5cblx0d2luZG93LnNwbi5vbmF1ZGlvcHJvY2VzcyA9IGZ1bmN0aW9uKGUpIHtcblx0XHQvLyBcInNlbGZcIiBsYW5kIG92ZXIgaGVyZSBiZWNhdXNlIG9mIHRoZSBmdW5jdGlvbiBjbG9zdXJlXG5cblx0XHQvL3RoaXMgaXMgdG8gb2J0YWluIHRoZSBjdXJyZW50IGFtcGxpdHVkZSBzcGVjdHJ1bVxuXHRcdHZhciBzaWduYWwgPSBzZWxmLnNpZ25hbCA9IGUuaW5wdXRCdWZmZXIuZ2V0Q2hhbm5lbERhdGEoMCk7XG5cdFx0dmFyIGRhdGEgPSBzZWxmLnNwZWN0cnVtRGF0YTtcblx0XHR2YXIgc3BlYyA9IHNlbGYuY29tcGxleFNwZWN0cnVtO1xuXHRcdHZhciBhbXBTcGVjdHJ1bSA9IHNlbGYuYW1wU3BlY3RydW07XG5cdFx0dmFyIHdpbmRvd2VkU2lnbmFsID0gc2VsZi5jb21wdXRlV2luZG93KHNpZ25hbCwgc2VsZi53aW5kb3dpbmdGdW5jdGlvbik7XG5cblx0XHQvL21hcCB0aW1lIGRvbWFpblxuXHRcdGRhdGEubWFwKGZ1bmN0aW9uKHZhbHVlLCBpLCBuKSB7XG5cdFx0XHR2YWx1ZS5yZWFsID0gd2luZG93ZWRTaWduYWxbaV07XG5cdFx0fSk7XG5cblx0XHQvL2NhbGN1bGF0ZSBhbXBsaXR1ZGVcblx0XHRzZWxmLmNvbXB1dGVBbXBsaXR1ZGUoc3BlYywgYW1wU3BlY3RydW0sIGJ1ZmZlclNpemUpO1xuXG5cdFx0Ly9jYWxsIGNhbGxiYWNrIGlmIGFwcGxpY2FibGVcblx0XHRpZiAodHlwZW9mIGNhbGxiYWNrID09PSBcImZ1bmN0aW9uXCIgJiYgRVhUUkFDVElPTl9TVEFSVEVEKSB7XG5cdFx0XHRjYWxsYmFjayhzZWxmLmdldChzZWxmLl9mZWF0dXJlc1RvRXh0cmFjdCkpO1xuXHRcdH1cblxuXHR9O1xuXG5cdHdpbmRvdy5zcG4uY29ubmVjdChhdWRpb0NvbnRleHQuZGVzdGluYXRpb24pO1xuXHRzb3VyY2UuY29ubmVjdCh3aW5kb3cuc3BuLCAwLCAwKTtcblxuXHQvLyBjb25zdHJ1Y3RvcnMgcmV0dXJuIFwidGhpc1wiIGJ5IGRlZmF1bHRcbn1cblxuXG4vLyBDb21wdXRlIG1ldGhvZHNcbi8vIC0tLS0tLS0tLS0tLS0tLVxuXG5NZXlkYS5wcm90b3R5cGUuY29tcHV0ZUFtcGxpdHVkZSA9IGZ1bmN0aW9uKGNvbXBsZXhTcGVjdHJ1bSwgYW1wU3BlY3RydW0sIGJ1ZmZlclNpemUpIHtcblxuXHQvLyB3b3JrcyBpbnNpZGUgYW5kIG91dHNpZGUgTWV5ZGFcblx0Y29tcGxleFNwZWN0cnVtID0gY29tcGxleFNwZWN0cnVtIHx8wqB0aGlzLmNvbXBsZXhTcGVjdHJ1bTtcblx0YW1wU3BlY3RydW0gXHRcdD0gYW1wU3BlY3RydW0gXHRcdHx8wqB0aGlzLmFtcFNwZWN0cnVtO1xuXHRidWZmZXJTaXplIFx0XHRcdD0gYnVmZmVyU2l6ZSBcdFx0XHR8fMKgdGhpcy5idWZmZXJTaXplO1xuXG5cdGZvciAodmFyIGkgPSAwOyBpIDwgYnVmZmVyU2l6ZS8yOyBpKyspIHtcblx0XHRhbXBTcGVjdHJ1bVtpXSA9IE1hdGguc3FydChNYXRoLnBvdyhjb21wbGV4U3BlY3RydW0ucmVhbFtpXSwyKSArIE1hdGgucG93KGNvbXBsZXhTcGVjdHJ1bS5pbWFnW2ldLDIpKTtcblx0fVxufTtcblxuTWV5ZGEucHJvdG90eXBlLmNvbXB1dGVIYW1taW5nID0gZnVuY3Rpb24oYnVmZmVyU2l6ZSkge1xuXHRidWZmZXJTaXplID0gYnVmZmVyU2l6ZSB8fMKgdGhpcy5idWZmZXJTaXplO1xuXG5cdHZhciBoYW1taW5nID0gbmV3IEZsb2F0MzJBcnJheShidWZmZXJTaXplKTtcblx0Zm9yICh2YXIgaSA9IDA7IGkgPCBidWZmZXJTaXplOyBpKyspIHtcblx0XHQvL0FjY29yZGluZyB0byBodHRwOi8vdWsubWF0aHdvcmtzLmNvbS9oZWxwL3NpZ25hbC9yZWYvaGFtbWluZy5odG1sXG5cdFx0aGFtbWluZ1tpXSA9IDAuNTQgLSAwLjQ2Kk1hdGguY29zKDIqTWF0aC5QSSooaS9idWZmZXJTaXplLTEpKTtcblx0fVxuXG5cdHJldHVybiBoYW1taW5nO1xufTtcblxuTWV5ZGEucHJvdG90eXBlLmNvbXB1dGVIYW5uaW5nID0gZnVuY3Rpb24oYnVmZmVyU2l6ZSkge1xuXHRidWZmZXJTaXplID0gYnVmZmVyU2l6ZSB8fMKgdGhpcy5idWZmZXJTaXplO1xuXG5cdHZhciBoYW5uaW5nID0gbmV3IEZsb2F0MzJBcnJheShidWZmZXJTaXplKTtcblx0Zm9yICh2YXIgaSA9IDA7IGkgPCBidWZmZXJTaXplOyBpKyspIHtcblx0XHQvL0FjY29yZGluZyB0byB0aGUgUiBkb2N1bWVudGF0aW9uIGh0dHA6Ly9yZ20ub2dhbGFiLm5ldC9SR00vUl9yZGZpbGU/Zj1HRU5FQXJlYWQvbWFuL2hhbm5pbmcud2luZG93LlJkJmQ9Ul9DQ1xuXHRcdGhhbm5pbmdbaV0gPSAwLjUgLSAwLjUqTWF0aC5jb3MoMipNYXRoLlBJKmkvKGJ1ZmZlclNpemUtMSkpO1xuXHR9XG5cblx0cmV0dXJuIGhhbm5pbmc7XG59O1xuXG4vL1VORklOSVNIRUQgLSBibGFja21hbiB3aW5kb3cgaW1wbGVtZW50YXRpb25cbi8qXG5NZXlkYS5wcm90b3R5cGUuY29tcHV0ZUJsYWNrbWFuID0gZnVuY3Rpb24oYnVmZmVyU2l6ZSkge1xuXHRidWZmZXJTaXplID0gYnVmZmVyU2l6ZSB8fMKgdGhpcy5idWZmZXJTaXplO1xuXHRcblx0dmFyIGJsYWNrbWFuID0gbmV3IEZsb2F0MzJBcnJheShidWZmZXJTaXplKTtcblx0Ly9BY2NvcmRpbmcgdG8gaHR0cDovL3VrLm1hdGh3b3Jrcy5jb20vaGVscC9zaWduYWwvcmVmL2JsYWNrbWFuLmh0bWxcblx0Ly9maXJzdCBoYWxmIG9mIHRoZSB3aW5kb3dcblx0Zm9yICh2YXIgaSA9IDA7IGkgPCAoYnVmZmVyU2l6ZSAlIDIpID8gKGJ1ZmZlclNpemUrMSkvMiA6IGJ1ZmZlclNpemUvMjsgaSsrKSB7XG5cdFx0dGhpcy5ibGFja21hbltpXSA9IDAuNDIgLSAwLjUqTWF0aC5jb3MoMipNYXRoLlBJKmkvKGJ1ZmZlclNpemUtMSkpICsgMC4wOCpNYXRoLmNvcyg0Kk1hdGguUEkqaS8oYnVmZmVyU2l6ZS0xKSk7XG5cdH1cblx0Ly9zZWNvbmQgaGFsZiBvZiB0aGUgd2luZG93XG5cdGZvciAodmFyIGkgPSBidWZmZXJTaXplLzI7IGkgPiAwOyBpLS0pIHtcblx0XHR0aGlzLmJsYWNrbWFuW2J1ZmZlclNpemUgLSBpXSA9IHRoaXMuYmxhY2ttYW5baV07XG5cdH1cbn07XG4qL1xuXG5NZXlkYS5wcm90b3R5cGUuY29tcHV0ZVdpbmRvdyA9IGZ1bmN0aW9uKHNpZywgdHlwZSkge1xuXG5cdHZhciBpLCBsZW4gPSBzaWcubGVuZ3RoO1xuXHR2YXIgd2luZG93ZWQgPSBuZXcgRmxvYXQzMkFycmF5KGxlbik7XG5cblx0Zm9yIChpID0gMDsgaSA8IGxlbjsgaSsrKSB7XG5cdFx0d2luZG93ZWRbaV0gPSBzaWdbaV0gKiB0aGlzW3R5cGVdW2ldO1xuXHR9XG5cdFxuXHRyZXR1cm4gd2luZG93ZWQ7XG59O1xuXG5NZXlkYS5wcm90b3R5cGUuY29tcHV0ZUJhcmtTY2FsZSA9IGZ1bmN0aW9uKGJ1ZmZlclNpemUsIHNhbXBsZVJhdGUpIHtcblx0YnVmZmVyU2l6ZSA9IGJ1ZmZlclNpemUgfHzCoHRoaXMuYnVmZmVyU2l6ZTtcblx0c2FtcGxlUmF0ZSA9IHNhbXBsZVJhdGUgfHzCoHRoaXMuc2FtcGxlUmF0ZTtcblxuICB2YXIgYmFya1NjYWxlID0gbmV3IEZsb2F0MzJBcnJheShidWZmZXJTaXplKTtcblxuICBmb3IodmFyIGkgPSAwOyBpIDwgYnVmZmVyU2l6ZTsgaSsrKXtcbiAgICBiYXJrU2NhbGVbaV0gPSBpICogc2FtcGxlUmF0ZSAvIChidWZmZXJTaXplKTtcbiAgICBiYXJrU2NhbGVbaV0gPSAxMyAqIE1hdGguYXRhbihiYXJrU2NhbGVbaV0vMTMxNS44KSArIDMuNSAqIE1hdGguYXRhbihNYXRoLnBvdygoYmFya1NjYWxlW2ldLzc1MTgpLCAyKSk7XG4gIH1cblxuICByZXR1cm4gYmFya1NjYWxlO1xufTtcblxuTWV5ZGEucHJvdG90eXBlLmNvbXB1dGVTcGVjdHJ1bURhdGEgPSBmdW5jdGlvbihidWZmZXJTaXplKSB7XG5cdGJ1ZmZlclNpemUgPSBidWZmZXJTaXplIHx8wqB0aGlzLmJ1ZmZlclNpemU7XG5cblx0Ly9jcmVhdGUgY29tcGxleGFycmF5IHRvIGhvbGQgdGhlIHNwZWN0cnVtXG5cdHZhciBkYXRhID0gbmV3IENvbXBsZXhBcnJheShidWZmZXJTaXplKTtcblx0dmFyIHNwZWN0cnVtID0gZGF0YS5GRlQoKTsgLy90cmFuc2Zvcm1cblx0dmFyIGFtcFNwZWN0cnVtID0gbmV3IEZsb2F0MzJBcnJheShidWZmZXJTaXplLzIpO1xuXG5cdHJldHVybiB7XG5cdFx0ZGF0YTogZGF0YSxcblx0XHRzcGVjdHJ1bTogc3BlY3RydW0sXG5cdFx0YW1wU3BlY3RydW06IGFtcFNwZWN0cnVtXG5cdH07XG59O1xuXG5cbi8vIE1leWRhIG1ldGhvZHNcbi8vIC0tLS0tLS0tLS0tLS1cblxuLy8gbG9hZHMgYWxsIHRoZSBleHRyYWN0b3Igb2JqZWN0c1xuLy8gaW5pdGlhbGl6ZXMgdGhlbSBhbmQgYmluZHMgdGhlbVxuLy8gdG8gdGhlIGZlYXR1cmVFeHRyYWN0b3JzXG4vLyBhbmQgZmVhdHVyZUluZm8gbGlzdHNcbi8vIEBOT1RFIChjL3NoKW91bGQgYmUgaGFuZGVsZWQgZGlmZmVyZW50bHlcbk1leWRhLnByb3RvdHlwZS5pbml0aWFsaXNlRXh0cmFjdG9ycyA9IGZ1bmN0aW9uKCkge1xuXG5cdHZhciBleHRyYWN0b3JzID0gcmVxdWlyZSgnLi9leHRyYWN0b3JzJyk7XG5cblx0Ly8gTG91ZG5lc3Ncblx0dmFyIGxvdWRuZXNzID0gZXh0cmFjdG9ycy5sb3VkbmVzcyh7XG5cdFx0TlVNX0JBUktfQkFORFM6IDI0LFxuXHRcdGJhcmtTY2FsZTogdGhpcy5iYXJrU2NhbGUsXG5cdFx0bm9ybWFsaXNlZFNwZWN0cnVtOiB0aGlzLmFtcFNwZWN0cnVtLFxuXHRcdHNhbXBsZVJhdGU6IHRoaXMuYXVkaW9Db250ZXh0LnNhbXBsZVJhdGVcblx0fSk7XG5cblx0dGhpcy5mZWF0dXJlRXh0cmFjdG9ycy5sb3VkbmVzcyA9IGxvdWRuZXNzO1xuXHR0aGlzLmZlYXR1cmVJbmZvLmxvdWRuZXNzID0gbG91ZG5lc3MuaW5mbztcblxuXHQvLyBSZXN0IG9mIHRoZSBleHRyYWN0b3JzXG5cdC8vIC4uLlxufTtcblxuLy9zb3VyY2Ugc2V0dGVyIG1ldGhvZFxuTWV5ZGEucHJvdG90eXBlLnNldFNvdXJjZSA9IGZ1bmN0aW9uKF9zcmMpIHtcblx0X3NyYy5jb25uZWN0KHdpbmRvdy5zcG4pO1xufTtcblxuTWV5ZGEucHJvdG90eXBlLnN0YXJ0ID0gZnVuY3Rpb24oZmVhdHVyZXMpIHtcblx0dGhpcy5fZmVhdHVyZXNUb0V4dHJhY3QgPSBmZWF0dXJlcztcblx0dGhpcy5FWFRSQUNUSU9OX1NUQVJURUQgPSB0cnVlO1xufTtcblxuTWV5ZGEucHJvdG90eXBlLnN0b3AgPSBmdW5jdGlvbigpIHtcblx0dGhpcy5fZmVhdHVyZXNUb0V4dHJhY3QgPSBudWxsO1xuXHR0aGlzLkVYVFJBQ1RJT05fU1RBUlRFRCA9IGZhbHNlO1xufTtcblxuLy8gZGF0YSBwdWxsaW5nXG5NZXlkYS5wcm90b3R5cGUuZ2V0ID0gZnVuY3Rpb24oZmVhdHVyZSkge1xuXG5cdGlmKHR5cGVvZiBmZWF0dXJlID09PSBcIm9iamVjdFwiKXtcblx0XHR2YXIgcmVzdWx0cyA9IHt9O1xuXHRcdGZvciAodmFyIHggPSAwOyB4IDwgZmVhdHVyZS5sZW5ndGg7IHgrKyl7XG5cdFx0XHR0cnl7XG5cdFx0XHRcdHJlc3VsdHNbZmVhdHVyZVt4XV0gPSAodGhpcy5mZWF0dXJlRXh0cmFjdG9yc1tmZWF0dXJlW3hdXS5wcm9jZXNzKHRoaXMuc2lnbmFsKSk7XG5cdFx0XHR9IGNhdGNoIChlKXtcblx0XHRcdFx0Y29uc29sZS5lcnJvcihlKTtcblx0XHRcdH1cblx0XHR9XG5cdFx0cmV0dXJuIHJlc3VsdHM7XG5cdH0gZWxzZSBpZiAodHlwZW9mIGZlYXR1cmUgPT09IFwic3RyaW5nXCIpe1xuXHRcdHJldHVybiB0aGlzLmZlYXR1cmVFeHRyYWN0b3JzW2ZlYXR1cmVdLnByb2Nlc3ModGhpcy5zaWduYWwpO1xuXHR9IGVsc2V7XG5cdFx0dGhyb3cgbmV3IEVycm9yKFwiSW52YWxpZCBGZWF0dXJlIEZvcm1hdFwiKTtcblx0fVxufTtcblxubW9kdWxlLmV4cG9ydHMgPSBNZXlkYTtcbiIsIm1vZHVsZS5leHBvcnRzLsK1ID0gZnVuY3Rpb24oaSwgYW1wbGl0dWRlU3BlY3Qpe1xuICB2YXIgbnVtZXJhdG9yID0gMDtcbiAgdmFyIGRlbm9taW5hdG9yID0gMDtcbiAgXG4gIGZvcih2YXIgayA9IDA7IGsgPCBhbXBsaXR1ZGVTcGVjdC5sZW5ndGg7IGsrKyl7XG4gICAgbnVtZXJhdG9yICs9IE1hdGgucG93KGssaSkqTWF0aC5hYnMoYW1wbGl0dWRlU3BlY3Rba10pO1xuICAgIGRlbm9taW5hdG9yICs9IGFtcGxpdHVkZVNwZWN0W2tdO1xuICB9XG5cbiAgcmV0dXJuIG51bWVyYXRvci9kZW5vbWluYXRvcjtcbn07XG5cbm1vZHVsZS5leHBvcnRzLmlzUG93ZXJPZlR3byA9IGZ1bmN0aW9uKG51bSkge1xuICB3aGlsZSAoKChudW0gJSAyKSA9PT0gMCkgJiYgbnVtID4gMSkge1xuICAgIG51bSAvPSAyO1xuICB9XG5cbiAgcmV0dXJuIChudW0gPT0gMSk7XG59OyJdfQ==

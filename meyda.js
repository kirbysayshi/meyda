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
  
  "buffer": function(bufferSize, m) {
    return m.signal;
  },

  "rms": function(bufferSize, m) {

    var rms = 0;
    for (var i = 0; i < m.signal.length; i++) {
      rms += Math.pow(m.signal[i], 2);
    }
    rms = rms / m.signal.length;
    rms = Math.sqrt(rms);

    return rms;
  },

  "energy": function(bufferSize, m) {
    var energy = 0;
    for (var i = 0; i < m.signal.length; i++) {
      energy += Math.pow(Math.abs(m.signal[i]), 2);
    }
    return energy;
  },

  "complexSpectrum": function(bufferSize, m) {
    return m.complexSpectrum;
  },

  "spectralSlope": function(bufferSize, m) {
    //linear regression
    var ampSum = 0;
    var freqSum = 0;
    var freqs = new Float32Array(m.ampSpectrum.length);
    var powFreqSum = 0;
    var ampFreqSum = 0;

    for (var i = 0; i < m.ampSpectrum.length; i++) {
      ampSum += m.ampSpectrum[i];
      var curFreq = i * m.audioContext.sampleRate / bufferSize;
      freqs[i] = curFreq;
      powFreqSum += curFreq * curFreq;
      freqSum += curFreq;
      ampFreqSum += curFreq * m.ampSpectrum[i];
    }
    return (m.ampSpectrum.length * ampFreqSum - freqSum * ampSum) / (ampSum * (powFreqSum - Math.pow(freqSum, 2)));
  },
 
  "spectralCentroid": function(bufferSize, m) {
    return µ(1, m.ampSpectrum);
  },
 
  "spectralRolloff": function(bufferSize, m) {
    var ampspec = m.ampSpectrum;
    //calculate nyquist bin
    var nyqBin = m.audioContext.sampleRate / (2 * (ampspec.length - 1));
    var ec = 0;
    for (var i = 0; i < ampspec.length; i++) {
      ec += ampspec[i];
    }
    var threshold = 0.99 * ec;
    var n = ampspec.length - 1;
    while (ec > threshold && n >= 0) {
      ec -= ampspec[n];
      --n;
    }
    return (n + 1) * nyqBin;
  },

  "spectralFlatness": function(bufferSize, m) {
    var ampspec = m.ampSpectrum;
    var numerator = 0;
    var denominator = 0;
    for (var i = 0; i < ampspec.length; i++) {
      numerator += Math.log(ampspec[i]);
      denominator += ampspec[i];
    }
    return Math.exp(numerator / ampspec.length) * ampspec.length / denominator;
  },

  "spectralSpread": function(bufferSize, m) {
    var ampspec = m.ampSpectrum;
    return Math.sqrt(µ(2, ampspec) - Math.pow(µ(1, ampspec), 2));
  },

  "spectralSkewness": function(bufferSize, m, spectrum) {
    var ampspec = m.ampSpectrum;
    var µ1 = µ(1, ampspec);
    var µ2 = µ(2, ampspec);
    var µ3 = µ(3, ampspec);
    var numerator = 2 * Math.pow(µ1, 3) - 3 * µ1 * µ2 + µ3;
    var denominator = Math.pow(Math.sqrt(µ2 - Math.pow(µ1, 2)), 3);
    return numerator / denominator;
  },

  "spectralKurtosis": function(bufferSize, m) {
    var ampspec = m.ampSpectrum;
    var µ1 = µ(1, ampspec);
    var µ2 = µ(2, ampspec);
    var µ3 = µ(3, ampspec);
    var µ4 = µ(4, ampspec);
    var numerator = -3 * Math.pow(µ1, 4) + 6 * µ1 * µ2 - 4 * µ1 * µ3 + µ4;
    var denominator = Math.pow(Math.sqrt(µ2 - Math.pow(µ1, 2)), 4);
    return numerator / denominator;
  },

  "amplitudeSpectrum": function(bufferSize, m) {
    return m.ampSpectrum;
  },

  "zcr": function(bufferSize, m) {
    var zcr = 0;
    for (var i = 0; i < m.signal.length; i++) {
      if ((m.signal[i] >= 0 && m.signal[i + 1] < 0) || (m.signal[i] < 0 && m.signal[i + 1] >= 0)) {
        zcr++;
      }
    }
    return zcr;
  },
  
  "powerSpectrum": function(bufferSize, m) {
    var powerSpectrum = new Float32Array(m.ampSpectrum.length);
    for (var i = 0; i < powerSpectrum.length; i++) {
      powerSpectrum[i] = Math.pow(m.ampSpectrum[i], 2);
    }
    return powerSpectrum;
  },
  
  "loudness": function(bufferSize, m) {

    var NUM_BARK_BANDS = 24;
    var specific = new Float32Array(NUM_BARK_BANDS);
    var tot = 0;
    var normalisedSpectrum = m.ampSpectrum;
    var bbLimits = new Int32Array(NUM_BARK_BANDS + 1);

    bbLimits[0] = 0;
    var currentBandEnd = self.barkScale[m.ampSpectrum.length - 1] / NUM_BARK_BANDS;
    var currentBand = 1;
    for (var i = 0; i < m.ampSpectrum.length; i++) {
      while (self.barkScale[i] > currentBandEnd) {
        bbLimits[currentBand++] = i;
        currentBandEnd = currentBand * self.barkScale[m.ampSpectrum.length - 1] / NUM_BARK_BANDS;
      }
    }

    bbLimits[NUM_BARK_BANDS] = m.ampSpectrum.length - 1;

    //process

    for (var i = 0; i < NUM_BARK_BANDS; i++) {
      var sum = 0;
      for (var j = bbLimits[i]; j < bbLimits[i + 1]; j++) {

        sum += normalisedSpectrum[j];
      }
      specific[i] = Math.pow(sum, 0.23);
    }

    //get total loudness
    for (var i = 0; i < specific.length; i++) {
      tot += specific[i];
    }
    return {
      "specific": specific,
      "total": tot
    };
  },
  
  "perceptualSpread": function(bufferSize, m) {
    var loudness = m.featureExtractors.loudness(bufferSize, m);

    var max = 0;
    for (var i = 0; i < loudness.specific.length; i++) {
      if (loudness.specific[i] > max) {
        max = loudness.specific[i];
      }
    }

    var spread = Math.pow((loudness.total - max) / loudness.total, 2);

    return spread;
  },
  
  "perceptualSharpness": function(bufferSize, m) {
    var loudness = m.featureExtractors.loudness(bufferSize, m);
    var spec = loudness.specific;
    var output = 0;

    for (var i = 0; i < spec.length; i++) {
      if (i < 15) {
        output += (i + 1) * spec[i + 1];
      } else {
        output += 0.066 * Math.exp(0.171 * (i + 1));
      }
    }
    output *= 0.11 / loudness.total;

    return output;
  },
  
  "mfcc": function(bufferSize, m) {
    //used tutorial from http://practicalcryptography.com/miscellaneous/machine-learning/guide-mel-frequency-cepstral-coefficients-mfccs/
    var powSpec = m.featureExtractors.powerSpectrum(bufferSize, m);
    var freqToMel = function(freqValue) {
      var melValue = 1125 * Math.log(1 + (freqValue / 700));
      return melValue;
    };
    var melToFreq = function(melValue) {
      var freqValue = 700 * (Math.exp(melValue / 1125) - 1);
      return freqValue;
    };
    var numFilters = 26; //26 filters is standard
    var melValues = new Float32Array(numFilters + 2); //the +2 is the upper and lower limits
    var melValuesInFreq = new Float32Array(numFilters + 2);
    //Generate limits in Hz - from 0 to the nyquist.
    var lowerLimitFreq = 0;
    var upperLimitFreq = audioContext.sampleRate / 2;
    //Convert the limits to Mel
    var lowerLimitMel = freqToMel(lowerLimitFreq);
    var upperLimitMel = freqToMel(upperLimitFreq);
    //Find the range
    var range = upperLimitMel - lowerLimitMel;
    //Find the range as part of the linear interpolation
    var valueToAdd = range / (numFilters + 1);

    var fftBinsOfFreq = Array(numFilters + 2);

    for (var i = 0; i < melValues.length; i++) {
      //Initialising the mel frequencies - they are just a linear interpolation between the lower and upper limits.
      melValues[i] = i * valueToAdd;
      //Convert back to Hz
      melValuesInFreq[i] = melToFreq(melValues[i]);
      //Find the corresponding bins
      fftBinsOfFreq[i] = Math.floor((bufferSize + 1) * melValuesInFreq[i] / audioContext.sampleRate);
    }

    var filterBank = Array(numFilters);
    for (var j = 0; j < filterBank.length; j++) {
      //creating a two dimensional array of size numFiltes * (buffersize/2)+1 and pre-populating the arrays with 0s.
      filterBank[j] = Array.apply(null, new Array((bufferSize / 2) + 1)).map(Number.prototype.valueOf, 0);
      //creating the lower and upper slopes for each bin
      for (var i = fftBinsOfFreq[j]; i < fftBinsOfFreq[j + 1]; i++) {
        filterBank[j][i] = (i - fftBinsOfFreq[j]) / (fftBinsOfFreq[j + 1] - fftBinsOfFreq[j]);
      }
      for (var i = fftBinsOfFreq[j + 1]; i < fftBinsOfFreq[j + 2]; i++) {
        filterBank[j][i] = (fftBinsOfFreq[j + 2] - i) / (fftBinsOfFreq[j + 2] - fftBinsOfFreq[j + 1]);
      }
    }

    var loggedMelBands = new Float32Array(numFilters);
    for (var i = 0; i < loggedMelBands.length; i++) {
      loggedMelBands[i] = 0;
      for (var j = 0; j < (bufferSize / 2); j++) {
        //point multiplication between power spectrum and filterbanks.
        filterBank[i][j] = filterBank[i][j] * powSpec[j];

        //summing up all of the coefficients into one array
        loggedMelBands[i] += filterBank[i][j];
      }
      //log each coefficient
      loggedMelBands[i] = Math.log(loggedMelBands[i]);
    }

    //dct
    var k = Math.PI / numFilters;
    var w1 = 1.0 / Math.sqrt(numFilters);
    var w2 = Math.sqrt(2.0 / numFilters);
    var numCoeffs = 13;
    var dctMatrix = new Float32Array(numCoeffs * numFilters);

    for (var i = 0; i < numCoeffs; i++) {
      for (var j = 0; j < numFilters; j++) {
        var idx = i + (j * numCoeffs);
        if (i == 0) {
          dctMatrix[idx] = w1 * Math.cos(k * (i + 1) * (j + 0.5));
        } else {
          dctMatrix[idx] = w2 * Math.cos(k * (i + 1) * (j + 0.5));
        }
      }
    }

    var mfccs = new Float32Array(numCoeffs);
    for (var k = 0; k < numCoeffs; k++) {
      var v = 0;
      for (var n = 0; n < numFilters; n++) {
        var idx = k + (n * numCoeffs);
        v += (dctMatrix[idx] * loggedMelBands[n]);
      }
      mfccs[k] = v / numCoeffs;
    }
    return mfccs;
  }
};
},{}],4:[function(require,module,exports){
// Meyda Javascript DSP library

var ComplexArray = require('../lib/jsfft/complex_array').ComplexArray;
// modifies ComplexArray
var fft = require('../lib/jsfft/fft');

module.exports = function(audioContext,src,bufSize,callback){
	
	//I am myself
	var self = this;
	
	self.featureExtractors = require('./extractors');

	//default buffer size
	var bufferSize = bufSize ? bufSize : 256;

	//initial source
	var source = src;

	//callback controllers
	var EXTRACTION_STARTED = false;
	var _featuresToExtract;

	//utilities
	var µ = function(i, amplitudeSpect){
		var numerator = 0;
		var denominator = 0;
		for(var k = 0; k < amplitudeSpect.length; k++){
			numerator += Math.pow(k,i)*Math.abs(amplitudeSpect[k]);
			denominator += amplitudeSpect[k];
		}
		return numerator/denominator;
	};

	var isPowerOfTwo = function(num) {
		while (((num % 2) == 0) && num > 1) {
			num /= 2;
		}
		return (num == 1);
	};

	//initilize bark scale (to be used in most perceptual features).
	self.barkScale = new Float32Array(bufSize);

	for(var i = 0; i < self.barkScale.length; i++){
		self.barkScale[i] = i*audioContext.sampleRate/(bufSize);
		self.barkScale[i] = 13*Math.atan(self.barkScale[i]/1315.8) + 3.5* Math.atan(Math.pow((self.barkScale[i]/7518),2));
	}

	//WINDOWING
	//set default
	self.windowingFunction = "hanning";

	//create windows
	self.hanning = new Float32Array(bufSize);
	for (var i = 0; i < bufSize; i++) {
		//According to the R documentation http://rgm.ogalab.net/RGM/R_rdfile?f=GENEAread/man/hanning.window.Rd&d=R_CC
		self.hanning[i] = 0.5 - 0.5*Math.cos(2*Math.PI*i/(bufSize-1));
	}

	self.hamming = new Float32Array(bufSize);
	for (var i = 0; i < bufSize; i++) {
		//According to http://uk.mathworks.com/help/signal/ref/hamming.html
		self.hamming[i] = 0.54 - 0.46*Math.cos(2*Math.PI*(i/bufSize-1));
	}

	//UNFINISHED - blackman window implementation

	/*self.blackman = new Float32Array(bufSize);
	//According to http://uk.mathworks.com/help/signal/ref/blackman.html
	//first half of the window
	for (var i = 0; i < (bufSize % 2) ? (bufSize+1)/2 : bufSize/2; i++) {
		self.blackman[i] = 0.42 - 0.5*Math.cos(2*Math.PI*i/(bufSize-1)) + 0.08*Math.cos(4*Math.PI*i/(bufSize-1));
	}
	//second half of the window
	for (var i = bufSize/2; i > 0; i--) {
		self.blackman[bufSize - i] = self.blackman[i];
	}*/

	self.windowing = function(sig, type){
		var windowed = new Float32Array(sig.length);
		var i ,len = sig.length;

		if (type == "hanning") {
			for (i = 0; i < len; i++) {
				windowed[i] = sig[i]*self.hanning[i];
			}
		}
		else if (type == "hamming") {
			for (i = 0; i < len; i++) {
				windowed[i] = sig[i]*self.hamming[i];
			}
		}
		else if (type == "blackman") {
			for (i = 0; i < len; i++) {
				windowed[i] = sig[i]*self.blackman[i];
			}
		}

		return windowed;
	};

	//source setter method
	self.setSource = function(_src) {
		source = _src;
		source.connect(window.spn);
	};



	if (isPowerOfTwo(bufferSize) && audioContext) {
			self.featureInfo = {
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

			//create complexarray to hold the spectrum
			var data = new ComplexArray(bufferSize);

			//create nodes
			window.spn = audioContext.createScriptProcessor(bufferSize,1,1);
			spn.connect(audioContext.destination);

			window.spn.onaudioprocess = function(e) {
				//this is to obtain the current amplitude spectrum
				var inputData = e.inputBuffer.getChannelData(0);
				self.signal = inputData;
				var windowedSignal = self.windowing(self.signal, self.windowingFunction);

				//map time domain
				data.map(function(value, i, n) {
					value.real = windowedSignal[i];
				});
				//transform
				var spec = data.FFT();
				//assign to meyda
				self.complexSpectrum = spec;
				self.ampSpectrum = new Float32Array(bufferSize/2);
				//calculate amplitude
				for (var i = 0; i < bufferSize/2; i++) {
					self.ampSpectrum[i] = Math.sqrt(Math.pow(spec.real[i],2) + Math.pow(spec.imag[i],2));

				}
				//call callback if applicable
				if (typeof callback === "function" && EXTRACTION_STARTED) {
					callback(self.get(_featuresToExtract));
				}

			};

			self.start = function(features) {
				_featuresToExtract = features;
				EXTRACTION_STARTED = true;
			};

			self.stop = function() {
				EXTRACTION_STARTED = false;
			};

			self.audioContext = audioContext;

			self.get = function(feature) {
				if(typeof feature === "object"){
					var results = {};
					for (var x = 0; x < feature.length; x++){
						try{
							results[feature[x]] = (self.featureExtractors[feature[x]](bufferSize, self));
						} catch (e){
							console.error(e);
						}
					}
					return results;
				}
				else if (typeof feature === "string"){
					return self.featureExtractors[feature](bufferSize, self);
				}
				else{
					throw "Invalid Feature Format";
				}
			};
			source.connect(window.spn, 0, 0);
			return self;
	}
	else {
		//handle errors
		if (typeof audioContext == "undefined") {
			throw "AudioContext wasn't specified: Meyda will not run.";
		}
		else {
			throw "Buffer size is not a power of two: Meyda will not run.";
		}
	}
};
},{"../lib/jsfft/complex_array":1,"../lib/jsfft/fft":2,"./extractors":3}]},{},[4])(4)
});
//# sourceMappingURL=data:application/json;base64,eyJ2ZXJzaW9uIjozLCJzb3VyY2VzIjpbIi4uLy4uLy4uLy4uL3Vzci9sb2NhbC9saWIvbm9kZV9tb2R1bGVzL2Jyb3dzZXJpZnkvbm9kZV9tb2R1bGVzL2Jyb3dzZXItcGFjay9fcHJlbHVkZS5qcyIsImxpYi9qc2ZmdC9jb21wbGV4X2FycmF5LmpzIiwibGliL2pzZmZ0L2ZmdC5qcyIsInNyYy9leHRyYWN0b3JzL2luZGV4LmpzIiwic3JjL21leWRhLmpzIl0sIm5hbWVzIjpbXSwibWFwcGluZ3MiOiJBQUFBO0FDQUE7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBOztBQ3BIQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTs7QUNqT0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTs7QUN2U0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0EiLCJmaWxlIjoiZ2VuZXJhdGVkLmpzIiwic291cmNlUm9vdCI6IiIsInNvdXJjZXNDb250ZW50IjpbIihmdW5jdGlvbiBlKHQsbixyKXtmdW5jdGlvbiBzKG8sdSl7aWYoIW5bb10pe2lmKCF0W29dKXt2YXIgYT10eXBlb2YgcmVxdWlyZT09XCJmdW5jdGlvblwiJiZyZXF1aXJlO2lmKCF1JiZhKXJldHVybiBhKG8sITApO2lmKGkpcmV0dXJuIGkobywhMCk7dmFyIGY9bmV3IEVycm9yKFwiQ2Fubm90IGZpbmQgbW9kdWxlICdcIitvK1wiJ1wiKTt0aHJvdyBmLmNvZGU9XCJNT0RVTEVfTk9UX0ZPVU5EXCIsZn12YXIgbD1uW29dPXtleHBvcnRzOnt9fTt0W29dWzBdLmNhbGwobC5leHBvcnRzLGZ1bmN0aW9uKGUpe3ZhciBuPXRbb11bMV1bZV07cmV0dXJuIHMobj9uOmUpfSxsLGwuZXhwb3J0cyxlLHQsbixyKX1yZXR1cm4gbltvXS5leHBvcnRzfXZhciBpPXR5cGVvZiByZXF1aXJlPT1cImZ1bmN0aW9uXCImJnJlcXVpcmU7Zm9yKHZhciBvPTA7bzxyLmxlbmd0aDtvKyspcyhyW29dKTtyZXR1cm4gc30pIiwiJ3VzZSBzdHJpY3QnO1xuXG4hZnVuY3Rpb24oZXhwb3J0cywgdW5kZWZpbmVkKSB7XG5cbiAgdmFyXG4gICAgLy8gSWYgdGhlIHR5cGVkIGFycmF5IGlzIHVuc3BlY2lmaWVkLCB1c2UgdGhpcy5cbiAgICBEZWZhdWx0QXJyYXlUeXBlID0gRmxvYXQzMkFycmF5LFxuICAgIC8vIFNpbXBsZSBtYXRoIGZ1bmN0aW9ucyB3ZSBuZWVkLlxuICAgIHNxcnQgPSBNYXRoLnNxcnQsXG4gICAgc3FyID0gZnVuY3Rpb24obnVtYmVyKSB7cmV0dXJuIE1hdGgucG93KG51bWJlciwgMil9LFxuICAgIC8vIEludGVybmFsIGNvbnZlbmllbmNlIGNvcGllcyBvZiB0aGUgZXhwb3J0ZWQgZnVuY3Rpb25zXG4gICAgaXNDb21wbGV4QXJyYXksXG4gICAgQ29tcGxleEFycmF5XG5cbiAgZXhwb3J0cy5pc0NvbXBsZXhBcnJheSA9IGlzQ29tcGxleEFycmF5ID0gZnVuY3Rpb24ob2JqKSB7XG4gICAgcmV0dXJuIG9iaiAhPT0gdW5kZWZpbmVkICYmXG4gICAgICBvYmouaGFzT3duUHJvcGVydHkgIT09IHVuZGVmaW5lZCAmJlxuICAgICAgb2JqLmhhc093blByb3BlcnR5KCdyZWFsJykgJiZcbiAgICAgIG9iai5oYXNPd25Qcm9wZXJ0eSgnaW1hZycpXG4gIH1cblxuICBleHBvcnRzLkNvbXBsZXhBcnJheSA9IENvbXBsZXhBcnJheSA9IGZ1bmN0aW9uKG90aGVyLCBvcHRfYXJyYXlfdHlwZSl7XG4gICAgaWYgKGlzQ29tcGxleEFycmF5KG90aGVyKSkge1xuICAgICAgLy8gQ29weSBjb25zdHVjdG9yLlxuICAgICAgdGhpcy5BcnJheVR5cGUgPSBvdGhlci5BcnJheVR5cGVcbiAgICAgIHRoaXMucmVhbCA9IG5ldyB0aGlzLkFycmF5VHlwZShvdGhlci5yZWFsKVxuICAgICAgdGhpcy5pbWFnID0gbmV3IHRoaXMuQXJyYXlUeXBlKG90aGVyLmltYWcpXG4gICAgfSBlbHNlIHtcbiAgICAgIHRoaXMuQXJyYXlUeXBlID0gb3B0X2FycmF5X3R5cGUgfHwgRGVmYXVsdEFycmF5VHlwZVxuICAgICAgLy8gb3RoZXIgY2FuIGJlIGVpdGhlciBhbiBhcnJheSBvciBhIG51bWJlci5cbiAgICAgIHRoaXMucmVhbCA9IG5ldyB0aGlzLkFycmF5VHlwZShvdGhlcilcbiAgICAgIHRoaXMuaW1hZyA9IG5ldyB0aGlzLkFycmF5VHlwZSh0aGlzLnJlYWwubGVuZ3RoKVxuICAgIH1cblxuICAgIHRoaXMubGVuZ3RoID0gdGhpcy5yZWFsLmxlbmd0aFxuICB9XG5cbiAgQ29tcGxleEFycmF5LnByb3RvdHlwZS50b1N0cmluZyA9IGZ1bmN0aW9uKCkge1xuICAgIHZhciBjb21wb25lbnRzID0gW11cblxuICAgIHRoaXMuZm9yRWFjaChmdW5jdGlvbihjX3ZhbHVlLCBpKSB7XG4gICAgICBjb21wb25lbnRzLnB1c2goXG4gICAgICAgICcoJyArXG4gICAgICAgIGNfdmFsdWUucmVhbC50b0ZpeGVkKDIpICsgJywnICtcbiAgICAgICAgY192YWx1ZS5pbWFnLnRvRml4ZWQoMikgK1xuICAgICAgICAnKSdcbiAgICAgIClcbiAgICB9KVxuXG4gICAgcmV0dXJuICdbJyArIGNvbXBvbmVudHMuam9pbignLCcpICsgJ10nXG4gIH1cblxuICAvLyBJbi1wbGFjZSBtYXBwZXIuXG4gIENvbXBsZXhBcnJheS5wcm90b3R5cGUubWFwID0gZnVuY3Rpb24obWFwcGVyKSB7XG4gICAgdmFyXG4gICAgICBpLFxuICAgICAgbiA9IHRoaXMubGVuZ3RoLFxuICAgICAgLy8gRm9yIEdDIGVmZmljaWVuY3ksIHBhc3MgYSBzaW5nbGUgY192YWx1ZSBvYmplY3QgdG8gdGhlIG1hcHBlci5cbiAgICAgIGNfdmFsdWUgPSB7fVxuXG4gICAgZm9yIChpID0gMDsgaSA8IG47IGkrKykge1xuICAgICAgY192YWx1ZS5yZWFsID0gdGhpcy5yZWFsW2ldXG4gICAgICBjX3ZhbHVlLmltYWcgPSB0aGlzLmltYWdbaV1cbiAgICAgIG1hcHBlcihjX3ZhbHVlLCBpLCBuKVxuICAgICAgdGhpcy5yZWFsW2ldID0gY192YWx1ZS5yZWFsXG4gICAgICB0aGlzLmltYWdbaV0gPSBjX3ZhbHVlLmltYWdcbiAgICB9XG5cbiAgICByZXR1cm4gdGhpc1xuICB9XG5cbiAgQ29tcGxleEFycmF5LnByb3RvdHlwZS5mb3JFYWNoID0gZnVuY3Rpb24oaXRlcmF0b3IpIHtcbiAgICB2YXJcbiAgICAgIGksXG4gICAgICBuID0gdGhpcy5sZW5ndGgsXG4gICAgICAvLyBGb3IgY29uc2lzdGVuY3kgd2l0aCAubWFwLlxuICAgICAgY192YWx1ZSA9IHt9XG5cbiAgICBmb3IgKGkgPSAwOyBpIDwgbjsgaSsrKSB7XG4gICAgICBjX3ZhbHVlLnJlYWwgPSB0aGlzLnJlYWxbaV1cbiAgICAgIGNfdmFsdWUuaW1hZyA9IHRoaXMuaW1hZ1tpXVxuICAgICAgaXRlcmF0b3IoY192YWx1ZSwgaSwgbilcbiAgICB9XG4gIH1cblxuICBDb21wbGV4QXJyYXkucHJvdG90eXBlLmNvbmp1Z2F0ZSA9IGZ1bmN0aW9uKCkge1xuICAgIHJldHVybiAobmV3IENvbXBsZXhBcnJheSh0aGlzKSkubWFwKGZ1bmN0aW9uKHZhbHVlKSB7XG4gICAgICB2YWx1ZS5pbWFnICo9IC0xXG4gICAgfSlcbiAgfVxuXG4gIC8vIEhlbHBlciBzbyB3ZSBjYW4gbWFrZSBBcnJheVR5cGUgb2JqZWN0cyByZXR1cm5lZCBoYXZlIHNpbWlsYXIgaW50ZXJmYWNlc1xuICAvLyAgIHRvIENvbXBsZXhBcnJheXMuXG4gIGZ1bmN0aW9uIGl0ZXJhYmxlKG9iaikge1xuICAgIGlmICghb2JqLmZvckVhY2gpXG4gICAgICBvYmouZm9yRWFjaCA9IGZ1bmN0aW9uKGl0ZXJhdG9yKSB7XG4gICAgICAgIHZhciBpLCBuID0gdGhpcy5sZW5ndGhcblxuICAgICAgICBmb3IgKGkgPSAwOyBpIDwgbjsgaSsrKVxuICAgICAgICAgIGl0ZXJhdG9yKHRoaXNbaV0sIGksIG4pXG4gICAgICB9XG5cbiAgICByZXR1cm4gb2JqXG4gIH1cblxuICBDb21wbGV4QXJyYXkucHJvdG90eXBlLm1hZ25pdHVkZSA9IGZ1bmN0aW9uKCkge1xuICAgIHZhciBtYWdzID0gbmV3IHRoaXMuQXJyYXlUeXBlKHRoaXMubGVuZ3RoKVxuXG4gICAgdGhpcy5mb3JFYWNoKGZ1bmN0aW9uKHZhbHVlLCBpKSB7XG4gICAgICBtYWdzW2ldID0gc3FydChzcXIodmFsdWUucmVhbCkgKyBzcXIodmFsdWUuaW1hZykpXG4gICAgfSlcblxuICAgIC8vIEFycmF5VHlwZSB3aWxsIG5vdCBuZWNlc3NhcmlseSBiZSBpdGVyYWJsZTogbWFrZSBpdCBzby5cbiAgICByZXR1cm4gaXRlcmFibGUobWFncylcbiAgfVxufSh0eXBlb2YgZXhwb3J0cyA9PT0gJ3VuZGVmaW5lZCcgJiYgKHRoaXMuY29tcGxleF9hcnJheSA9IHt9KSB8fCBleHBvcnRzKVxuIiwiJ3VzZSBzdHJpY3QnO1xuXG4hZnVuY3Rpb24oZXhwb3J0cywgY29tcGxleF9hcnJheSkge1xuXG4gIHZhclxuICAgIENvbXBsZXhBcnJheSA9IGNvbXBsZXhfYXJyYXkuQ29tcGxleEFycmF5LFxuICAgIC8vIE1hdGggY29uc3RhbnRzIGFuZCBmdW5jdGlvbnMgd2UgbmVlZC5cbiAgICBQSSA9IE1hdGguUEksXG4gICAgU1FSVDFfMiA9IE1hdGguU1FSVDFfMixcbiAgICBzcXJ0ID0gTWF0aC5zcXJ0LFxuICAgIGNvcyA9IE1hdGguY29zLFxuICAgIHNpbiA9IE1hdGguc2luXG5cbiAgQ29tcGxleEFycmF5LnByb3RvdHlwZS5GRlQgPSBmdW5jdGlvbigpIHtcbiAgICByZXR1cm4gRkZUKHRoaXMsIGZhbHNlKTtcbiAgfVxuXG4gIGV4cG9ydHMuRkZUID0gZnVuY3Rpb24oaW5wdXQpIHtcbiAgICByZXR1cm4gZW5zdXJlQ29tcGxleEFycmF5KGlucHV0KS5GRlQoKVxuICB9XG5cbiAgQ29tcGxleEFycmF5LnByb3RvdHlwZS5JbnZGRlQgPSBmdW5jdGlvbigpIHtcbiAgICByZXR1cm4gRkZUKHRoaXMsIHRydWUpXG4gIH1cblxuICBleHBvcnRzLkludkZGVCA9IGZ1bmN0aW9uKGlucHV0KSB7XG4gICAgcmV0dXJuIGVuc3VyZUNvbXBsZXhBcnJheShpbnB1dCkuSW52RkZUKClcbiAgfVxuXG4gIC8vIEFwcGxpZXMgYSBmcmVxdWVuY3ktc3BhY2UgZmlsdGVyIHRvIGlucHV0LCBhbmQgcmV0dXJucyB0aGUgcmVhbC1zcGFjZVxuICAvLyBmaWx0ZXJlZCBpbnB1dC5cbiAgLy8gZmlsdGVyZXIgYWNjZXB0cyBmcmVxLCBpLCBuIGFuZCBtb2RpZmllcyBmcmVxLnJlYWwgYW5kIGZyZXEuaW1hZy5cbiAgQ29tcGxleEFycmF5LnByb3RvdHlwZS5mcmVxdWVuY3lNYXAgPSBmdW5jdGlvbihmaWx0ZXJlcikge1xuICAgIHJldHVybiB0aGlzLkZGVCgpLm1hcChmaWx0ZXJlcikuSW52RkZUKClcbiAgfVxuXG4gIGV4cG9ydHMuZnJlcXVlbmN5TWFwID0gZnVuY3Rpb24oaW5wdXQsIGZpbHRlcmVyKSB7XG4gICAgcmV0dXJuIGVuc3VyZUNvbXBsZXhBcnJheShpbnB1dCkuZnJlcXVlbmN5TWFwKGZpbHRlcmVyKVxuICB9XG5cbiAgZnVuY3Rpb24gZW5zdXJlQ29tcGxleEFycmF5KGlucHV0KSB7XG4gICAgcmV0dXJuIGNvbXBsZXhfYXJyYXkuaXNDb21wbGV4QXJyYXkoaW5wdXQpICYmIGlucHV0IHx8XG4gICAgICAgIG5ldyBDb21wbGV4QXJyYXkoaW5wdXQpXG4gIH1cblxuICBmdW5jdGlvbiBGRlQoaW5wdXQsIGludmVyc2UpIHtcbiAgICB2YXIgbiA9IGlucHV0Lmxlbmd0aFxuXG4gICAgaWYgKG4gJiAobiAtIDEpKSB7XG4gICAgICByZXR1cm4gRkZUX1JlY3Vyc2l2ZShpbnB1dCwgaW52ZXJzZSlcbiAgICB9IGVsc2Uge1xuICAgICAgcmV0dXJuIEZGVF8yX0l0ZXJhdGl2ZShpbnB1dCwgaW52ZXJzZSlcbiAgICB9XG4gIH1cblxuICBmdW5jdGlvbiBGRlRfUmVjdXJzaXZlKGlucHV0LCBpbnZlcnNlKSB7XG4gICAgdmFyXG4gICAgICBuID0gaW5wdXQubGVuZ3RoLFxuICAgICAgLy8gQ291bnRlcnMuXG4gICAgICBpLCBqLFxuICAgICAgb3V0cHV0LFxuICAgICAgLy8gQ29tcGxleCBtdWx0aXBsaWVyIGFuZCBpdHMgZGVsdGEuXG4gICAgICBmX3IsIGZfaSwgZGVsX2ZfciwgZGVsX2ZfaSxcbiAgICAgIC8vIExvd2VzdCBkaXZpc29yIGFuZCByZW1haW5kZXIuXG4gICAgICBwLCBtLFxuICAgICAgbm9ybWFsaXNhdGlvbixcbiAgICAgIHJlY3Vyc2l2ZV9yZXN1bHQsXG4gICAgICBfc3dhcCwgX3JlYWwsIF9pbWFnXG5cbiAgICBpZiAobiA9PT0gMSkge1xuICAgICAgcmV0dXJuIGlucHV0XG4gICAgfVxuXG4gICAgb3V0cHV0ID0gbmV3IENvbXBsZXhBcnJheShuLCBpbnB1dC5BcnJheVR5cGUpXG5cbiAgICAvLyBVc2UgdGhlIGxvd2VzdCBvZGQgZmFjdG9yLCBzbyB3ZSBhcmUgYWJsZSB0byB1c2UgRkZUXzJfSXRlcmF0aXZlIGluIHRoZVxuICAgIC8vIHJlY3Vyc2l2ZSB0cmFuc2Zvcm1zIG9wdGltYWxseS5cbiAgICBwID0gTG93ZXN0T2RkRmFjdG9yKG4pXG4gICAgbSA9IG4gLyBwXG4gICAgbm9ybWFsaXNhdGlvbiA9IDEgLyBzcXJ0KHApXG4gICAgcmVjdXJzaXZlX3Jlc3VsdCA9IG5ldyBDb21wbGV4QXJyYXkobSwgaW5wdXQuQXJyYXlUeXBlKVxuXG4gICAgLy8gTG9vcHMgZ28gbGlrZSBPKG4gzqMgcF9pKSwgd2hlcmUgcF9pIGFyZSB0aGUgcHJpbWUgZmFjdG9ycyBvZiBuLlxuICAgIC8vIGZvciBhIHBvd2VyIG9mIGEgcHJpbWUsIHAsIHRoaXMgcmVkdWNlcyB0byBPKG4gcCBsb2dfcCBuKVxuICAgIGZvcihqID0gMDsgaiA8IHA7IGorKykge1xuICAgICAgZm9yKGkgPSAwOyBpIDwgbTsgaSsrKSB7XG4gICAgICAgIHJlY3Vyc2l2ZV9yZXN1bHQucmVhbFtpXSA9IGlucHV0LnJlYWxbaSAqIHAgKyBqXVxuICAgICAgICByZWN1cnNpdmVfcmVzdWx0LmltYWdbaV0gPSBpbnB1dC5pbWFnW2kgKiBwICsgal1cbiAgICAgIH1cbiAgICAgIC8vIERvbid0IGdvIGRlZXBlciB1bmxlc3MgbmVjZXNzYXJ5IHRvIHNhdmUgYWxsb2NzLlxuICAgICAgaWYgKG0gPiAxKSB7XG4gICAgICAgIHJlY3Vyc2l2ZV9yZXN1bHQgPSBGRlQocmVjdXJzaXZlX3Jlc3VsdCwgaW52ZXJzZSlcbiAgICAgIH1cblxuICAgICAgZGVsX2ZfciA9IGNvcygyKlBJKmovbilcbiAgICAgIGRlbF9mX2kgPSAoaW52ZXJzZSA/IC0xIDogMSkgKiBzaW4oMipQSSpqL24pXG4gICAgICBmX3IgPSAxXG4gICAgICBmX2kgPSAwXG5cbiAgICAgIGZvcihpID0gMDsgaSA8IG47IGkrKykge1xuICAgICAgICBfcmVhbCA9IHJlY3Vyc2l2ZV9yZXN1bHQucmVhbFtpICUgbV1cbiAgICAgICAgX2ltYWcgPSByZWN1cnNpdmVfcmVzdWx0LmltYWdbaSAlIG1dXG5cbiAgICAgICAgb3V0cHV0LnJlYWxbaV0gKz0gZl9yICogX3JlYWwgLSBmX2kgKiBfaW1hZ1xuICAgICAgICBvdXRwdXQuaW1hZ1tpXSArPSBmX3IgKiBfaW1hZyArIGZfaSAqIF9yZWFsXG5cbiAgICAgICAgX3N3YXAgPSBmX3IgKiBkZWxfZl9yIC0gZl9pICogZGVsX2ZfaVxuICAgICAgICBmX2kgPSBmX3IgKiBkZWxfZl9pICsgZl9pICogZGVsX2ZfclxuICAgICAgICBmX3IgPSBfc3dhcFxuICAgICAgfVxuICAgIH1cblxuICAgIC8vIENvcHkgYmFjayB0byBpbnB1dCB0byBtYXRjaCBGRlRfMl9JdGVyYXRpdmUgaW4tcGxhY2VuZXNzXG4gICAgLy8gVE9ETzogZmFzdGVyIHdheSBvZiBtYWtpbmcgdGhpcyBpbi1wbGFjZT9cbiAgICBmb3IoaSA9IDA7IGkgPCBuOyBpKyspIHtcbiAgICAgIGlucHV0LnJlYWxbaV0gPSBub3JtYWxpc2F0aW9uICogb3V0cHV0LnJlYWxbaV1cbiAgICAgIGlucHV0LmltYWdbaV0gPSBub3JtYWxpc2F0aW9uICogb3V0cHV0LmltYWdbaV1cbiAgICB9XG5cbiAgICByZXR1cm4gaW5wdXRcbiAgfVxuXG4gIGZ1bmN0aW9uIEZGVF8yX0l0ZXJhdGl2ZShpbnB1dCwgaW52ZXJzZSkge1xuICAgIHZhclxuICAgICAgbiA9IGlucHV0Lmxlbmd0aCxcbiAgICAgIC8vIENvdW50ZXJzLlxuICAgICAgaSwgaixcbiAgICAgIG91dHB1dCwgb3V0cHV0X3IsIG91dHB1dF9pLFxuICAgICAgLy8gQ29tcGxleCBtdWx0aXBsaWVyIGFuZCBpdHMgZGVsdGEuXG4gICAgICBmX3IsIGZfaSwgZGVsX2ZfciwgZGVsX2ZfaSwgdGVtcCxcbiAgICAgIC8vIFRlbXBvcmFyeSBsb29wIHZhcmlhYmxlcy5cbiAgICAgIGxfaW5kZXgsIHJfaW5kZXgsXG4gICAgICBsZWZ0X3IsIGxlZnRfaSwgcmlnaHRfciwgcmlnaHRfaSxcbiAgICAgIC8vIHdpZHRoIG9mIGVhY2ggc3ViLWFycmF5IGZvciB3aGljaCB3ZSdyZSBpdGVyYXRpdmVseSBjYWxjdWxhdGluZyBGRlQuXG4gICAgICB3aWR0aFxuXG4gICAgb3V0cHV0ID0gQml0UmV2ZXJzZUNvbXBsZXhBcnJheShpbnB1dClcbiAgICBvdXRwdXRfciA9IG91dHB1dC5yZWFsXG4gICAgb3V0cHV0X2kgPSBvdXRwdXQuaW1hZ1xuICAgIC8vIExvb3BzIGdvIGxpa2UgTyhuIGxvZyBuKTpcbiAgICAvLyAgIHdpZHRoIH4gbG9nIG47IGksaiB+IG5cbiAgICB3aWR0aCA9IDFcbiAgICB3aGlsZSAod2lkdGggPCBuKSB7XG4gICAgICBkZWxfZl9yID0gY29zKFBJL3dpZHRoKVxuICAgICAgZGVsX2ZfaSA9IChpbnZlcnNlID8gLTEgOiAxKSAqIHNpbihQSS93aWR0aClcbiAgICAgIGZvciAoaSA9IDA7IGkgPCBuLygyKndpZHRoKTsgaSsrKSB7XG4gICAgICAgIGZfciA9IDFcbiAgICAgICAgZl9pID0gMFxuICAgICAgICBmb3IgKGogPSAwOyBqIDwgd2lkdGg7IGorKykge1xuICAgICAgICAgIGxfaW5kZXggPSAyKmkqd2lkdGggKyBqXG4gICAgICAgICAgcl9pbmRleCA9IGxfaW5kZXggKyB3aWR0aFxuXG4gICAgICAgICAgbGVmdF9yID0gb3V0cHV0X3JbbF9pbmRleF1cbiAgICAgICAgICBsZWZ0X2kgPSBvdXRwdXRfaVtsX2luZGV4XVxuICAgICAgICAgIHJpZ2h0X3IgPSBmX3IgKiBvdXRwdXRfcltyX2luZGV4XSAtIGZfaSAqIG91dHB1dF9pW3JfaW5kZXhdXG4gICAgICAgICAgcmlnaHRfaSA9IGZfaSAqIG91dHB1dF9yW3JfaW5kZXhdICsgZl9yICogb3V0cHV0X2lbcl9pbmRleF1cblxuICAgICAgICAgIG91dHB1dF9yW2xfaW5kZXhdID0gU1FSVDFfMiAqIChsZWZ0X3IgKyByaWdodF9yKVxuICAgICAgICAgIG91dHB1dF9pW2xfaW5kZXhdID0gU1FSVDFfMiAqIChsZWZ0X2kgKyByaWdodF9pKVxuICAgICAgICAgIG91dHB1dF9yW3JfaW5kZXhdID0gU1FSVDFfMiAqIChsZWZ0X3IgLSByaWdodF9yKVxuICAgICAgICAgIG91dHB1dF9pW3JfaW5kZXhdID0gU1FSVDFfMiAqIChsZWZ0X2kgLSByaWdodF9pKVxuICAgICAgICAgIHRlbXAgPSBmX3IgKiBkZWxfZl9yIC0gZl9pICogZGVsX2ZfaVxuICAgICAgICAgIGZfaSA9IGZfciAqIGRlbF9mX2kgKyBmX2kgKiBkZWxfZl9yXG4gICAgICAgICAgZl9yID0gdGVtcFxuICAgICAgICB9XG4gICAgICB9XG4gICAgICB3aWR0aCA8PD0gMVxuICAgIH1cblxuICAgIHJldHVybiBvdXRwdXRcbiAgfVxuXG4gIGZ1bmN0aW9uIEJpdFJldmVyc2VJbmRleChpbmRleCwgbikge1xuICAgIHZhciBiaXRyZXZlcnNlZF9pbmRleCA9IDBcblxuICAgIHdoaWxlIChuID4gMSkge1xuICAgICAgYml0cmV2ZXJzZWRfaW5kZXggPDw9IDFcbiAgICAgIGJpdHJldmVyc2VkX2luZGV4ICs9IGluZGV4ICYgMVxuICAgICAgaW5kZXggPj49IDFcbiAgICAgIG4gPj49IDFcbiAgICB9XG4gICAgcmV0dXJuIGJpdHJldmVyc2VkX2luZGV4XG4gIH1cblxuICBmdW5jdGlvbiBCaXRSZXZlcnNlQ29tcGxleEFycmF5KGFycmF5KSB7XG4gICAgdmFyIG4gPSBhcnJheS5sZW5ndGgsXG4gICAgICAgIGZsaXBzID0ge30sXG4gICAgICAgIHN3YXAsXG4gICAgICAgIGlcblxuICAgIGZvcihpID0gMDsgaSA8IG47IGkrKykge1xuICAgICAgdmFyIHJfaSA9IEJpdFJldmVyc2VJbmRleChpLCBuKVxuXG4gICAgICBpZiAoZmxpcHMuaGFzT3duUHJvcGVydHkoaSkgfHwgZmxpcHMuaGFzT3duUHJvcGVydHkocl9pKSkgY29udGludWVcblxuICAgICAgc3dhcCA9IGFycmF5LnJlYWxbcl9pXVxuICAgICAgYXJyYXkucmVhbFtyX2ldID0gYXJyYXkucmVhbFtpXVxuICAgICAgYXJyYXkucmVhbFtpXSA9IHN3YXBcblxuICAgICAgc3dhcCA9IGFycmF5LmltYWdbcl9pXVxuICAgICAgYXJyYXkuaW1hZ1tyX2ldID0gYXJyYXkuaW1hZ1tpXVxuICAgICAgYXJyYXkuaW1hZ1tpXSA9IHN3YXBcblxuICAgICAgZmxpcHNbaV0gPSBmbGlwc1tyX2ldID0gdHJ1ZVxuICAgIH1cblxuICAgIHJldHVybiBhcnJheVxuICB9XG5cbiAgZnVuY3Rpb24gTG93ZXN0T2RkRmFjdG9yKG4pIHtcbiAgICB2YXIgZmFjdG9yID0gMyxcbiAgICAgICAgc3FydF9uID0gc3FydChuKVxuXG4gICAgd2hpbGUoZmFjdG9yIDw9IHNxcnRfbikge1xuICAgICAgaWYgKG4gJSBmYWN0b3IgPT09IDApIHJldHVybiBmYWN0b3JcbiAgICAgIGZhY3RvciA9IGZhY3RvciArIDJcbiAgICB9XG4gICAgcmV0dXJuIG5cbiAgfVxuXG59KFxuICB0eXBlb2YgZXhwb3J0cyA9PT0gJ3VuZGVmaW5lZCcgJiYgKHRoaXMuZmZ0ID0ge30pIHx8IGV4cG9ydHMsXG4gIHR5cGVvZiByZXF1aXJlID09PSAndW5kZWZpbmVkJyAmJiAodGhpcy5jb21wbGV4X2FycmF5KSB8fFxuICAgIHJlcXVpcmUoJy4vY29tcGxleF9hcnJheScpXG4pXG4iLCJcbm1vZHVsZS5leHBvcnRzID0ge1xuICBcbiAgXCJidWZmZXJcIjogZnVuY3Rpb24oYnVmZmVyU2l6ZSwgbSkge1xuICAgIHJldHVybiBtLnNpZ25hbDtcbiAgfSxcblxuICBcInJtc1wiOiBmdW5jdGlvbihidWZmZXJTaXplLCBtKSB7XG5cbiAgICB2YXIgcm1zID0gMDtcbiAgICBmb3IgKHZhciBpID0gMDsgaSA8IG0uc2lnbmFsLmxlbmd0aDsgaSsrKSB7XG4gICAgICBybXMgKz0gTWF0aC5wb3cobS5zaWduYWxbaV0sIDIpO1xuICAgIH1cbiAgICBybXMgPSBybXMgLyBtLnNpZ25hbC5sZW5ndGg7XG4gICAgcm1zID0gTWF0aC5zcXJ0KHJtcyk7XG5cbiAgICByZXR1cm4gcm1zO1xuICB9LFxuXG4gIFwiZW5lcmd5XCI6IGZ1bmN0aW9uKGJ1ZmZlclNpemUsIG0pIHtcbiAgICB2YXIgZW5lcmd5ID0gMDtcbiAgICBmb3IgKHZhciBpID0gMDsgaSA8IG0uc2lnbmFsLmxlbmd0aDsgaSsrKSB7XG4gICAgICBlbmVyZ3kgKz0gTWF0aC5wb3coTWF0aC5hYnMobS5zaWduYWxbaV0pLCAyKTtcbiAgICB9XG4gICAgcmV0dXJuIGVuZXJneTtcbiAgfSxcblxuICBcImNvbXBsZXhTcGVjdHJ1bVwiOiBmdW5jdGlvbihidWZmZXJTaXplLCBtKSB7XG4gICAgcmV0dXJuIG0uY29tcGxleFNwZWN0cnVtO1xuICB9LFxuXG4gIFwic3BlY3RyYWxTbG9wZVwiOiBmdW5jdGlvbihidWZmZXJTaXplLCBtKSB7XG4gICAgLy9saW5lYXIgcmVncmVzc2lvblxuICAgIHZhciBhbXBTdW0gPSAwO1xuICAgIHZhciBmcmVxU3VtID0gMDtcbiAgICB2YXIgZnJlcXMgPSBuZXcgRmxvYXQzMkFycmF5KG0uYW1wU3BlY3RydW0ubGVuZ3RoKTtcbiAgICB2YXIgcG93RnJlcVN1bSA9IDA7XG4gICAgdmFyIGFtcEZyZXFTdW0gPSAwO1xuXG4gICAgZm9yICh2YXIgaSA9IDA7IGkgPCBtLmFtcFNwZWN0cnVtLmxlbmd0aDsgaSsrKSB7XG4gICAgICBhbXBTdW0gKz0gbS5hbXBTcGVjdHJ1bVtpXTtcbiAgICAgIHZhciBjdXJGcmVxID0gaSAqIG0uYXVkaW9Db250ZXh0LnNhbXBsZVJhdGUgLyBidWZmZXJTaXplO1xuICAgICAgZnJlcXNbaV0gPSBjdXJGcmVxO1xuICAgICAgcG93RnJlcVN1bSArPSBjdXJGcmVxICogY3VyRnJlcTtcbiAgICAgIGZyZXFTdW0gKz0gY3VyRnJlcTtcbiAgICAgIGFtcEZyZXFTdW0gKz0gY3VyRnJlcSAqIG0uYW1wU3BlY3RydW1baV07XG4gICAgfVxuICAgIHJldHVybiAobS5hbXBTcGVjdHJ1bS5sZW5ndGggKiBhbXBGcmVxU3VtIC0gZnJlcVN1bSAqIGFtcFN1bSkgLyAoYW1wU3VtICogKHBvd0ZyZXFTdW0gLSBNYXRoLnBvdyhmcmVxU3VtLCAyKSkpO1xuICB9LFxuIFxuICBcInNwZWN0cmFsQ2VudHJvaWRcIjogZnVuY3Rpb24oYnVmZmVyU2l6ZSwgbSkge1xuICAgIHJldHVybiDCtSgxLCBtLmFtcFNwZWN0cnVtKTtcbiAgfSxcbiBcbiAgXCJzcGVjdHJhbFJvbGxvZmZcIjogZnVuY3Rpb24oYnVmZmVyU2l6ZSwgbSkge1xuICAgIHZhciBhbXBzcGVjID0gbS5hbXBTcGVjdHJ1bTtcbiAgICAvL2NhbGN1bGF0ZSBueXF1aXN0IGJpblxuICAgIHZhciBueXFCaW4gPSBtLmF1ZGlvQ29udGV4dC5zYW1wbGVSYXRlIC8gKDIgKiAoYW1wc3BlYy5sZW5ndGggLSAxKSk7XG4gICAgdmFyIGVjID0gMDtcbiAgICBmb3IgKHZhciBpID0gMDsgaSA8IGFtcHNwZWMubGVuZ3RoOyBpKyspIHtcbiAgICAgIGVjICs9IGFtcHNwZWNbaV07XG4gICAgfVxuICAgIHZhciB0aHJlc2hvbGQgPSAwLjk5ICogZWM7XG4gICAgdmFyIG4gPSBhbXBzcGVjLmxlbmd0aCAtIDE7XG4gICAgd2hpbGUgKGVjID4gdGhyZXNob2xkICYmIG4gPj0gMCkge1xuICAgICAgZWMgLT0gYW1wc3BlY1tuXTtcbiAgICAgIC0tbjtcbiAgICB9XG4gICAgcmV0dXJuIChuICsgMSkgKiBueXFCaW47XG4gIH0sXG5cbiAgXCJzcGVjdHJhbEZsYXRuZXNzXCI6IGZ1bmN0aW9uKGJ1ZmZlclNpemUsIG0pIHtcbiAgICB2YXIgYW1wc3BlYyA9IG0uYW1wU3BlY3RydW07XG4gICAgdmFyIG51bWVyYXRvciA9IDA7XG4gICAgdmFyIGRlbm9taW5hdG9yID0gMDtcbiAgICBmb3IgKHZhciBpID0gMDsgaSA8IGFtcHNwZWMubGVuZ3RoOyBpKyspIHtcbiAgICAgIG51bWVyYXRvciArPSBNYXRoLmxvZyhhbXBzcGVjW2ldKTtcbiAgICAgIGRlbm9taW5hdG9yICs9IGFtcHNwZWNbaV07XG4gICAgfVxuICAgIHJldHVybiBNYXRoLmV4cChudW1lcmF0b3IgLyBhbXBzcGVjLmxlbmd0aCkgKiBhbXBzcGVjLmxlbmd0aCAvIGRlbm9taW5hdG9yO1xuICB9LFxuXG4gIFwic3BlY3RyYWxTcHJlYWRcIjogZnVuY3Rpb24oYnVmZmVyU2l6ZSwgbSkge1xuICAgIHZhciBhbXBzcGVjID0gbS5hbXBTcGVjdHJ1bTtcbiAgICByZXR1cm4gTWF0aC5zcXJ0KMK1KDIsIGFtcHNwZWMpIC0gTWF0aC5wb3cowrUoMSwgYW1wc3BlYyksIDIpKTtcbiAgfSxcblxuICBcInNwZWN0cmFsU2tld25lc3NcIjogZnVuY3Rpb24oYnVmZmVyU2l6ZSwgbSwgc3BlY3RydW0pIHtcbiAgICB2YXIgYW1wc3BlYyA9IG0uYW1wU3BlY3RydW07XG4gICAgdmFyIMK1MSA9IMK1KDEsIGFtcHNwZWMpO1xuICAgIHZhciDCtTIgPSDCtSgyLCBhbXBzcGVjKTtcbiAgICB2YXIgwrUzID0gwrUoMywgYW1wc3BlYyk7XG4gICAgdmFyIG51bWVyYXRvciA9IDIgKiBNYXRoLnBvdyjCtTEsIDMpIC0gMyAqIMK1MSAqIMK1MiArIMK1MztcbiAgICB2YXIgZGVub21pbmF0b3IgPSBNYXRoLnBvdyhNYXRoLnNxcnQowrUyIC0gTWF0aC5wb3cowrUxLCAyKSksIDMpO1xuICAgIHJldHVybiBudW1lcmF0b3IgLyBkZW5vbWluYXRvcjtcbiAgfSxcblxuICBcInNwZWN0cmFsS3VydG9zaXNcIjogZnVuY3Rpb24oYnVmZmVyU2l6ZSwgbSkge1xuICAgIHZhciBhbXBzcGVjID0gbS5hbXBTcGVjdHJ1bTtcbiAgICB2YXIgwrUxID0gwrUoMSwgYW1wc3BlYyk7XG4gICAgdmFyIMK1MiA9IMK1KDIsIGFtcHNwZWMpO1xuICAgIHZhciDCtTMgPSDCtSgzLCBhbXBzcGVjKTtcbiAgICB2YXIgwrU0ID0gwrUoNCwgYW1wc3BlYyk7XG4gICAgdmFyIG51bWVyYXRvciA9IC0zICogTWF0aC5wb3cowrUxLCA0KSArIDYgKiDCtTEgKiDCtTIgLSA0ICogwrUxICogwrUzICsgwrU0O1xuICAgIHZhciBkZW5vbWluYXRvciA9IE1hdGgucG93KE1hdGguc3FydCjCtTIgLSBNYXRoLnBvdyjCtTEsIDIpKSwgNCk7XG4gICAgcmV0dXJuIG51bWVyYXRvciAvIGRlbm9taW5hdG9yO1xuICB9LFxuXG4gIFwiYW1wbGl0dWRlU3BlY3RydW1cIjogZnVuY3Rpb24oYnVmZmVyU2l6ZSwgbSkge1xuICAgIHJldHVybiBtLmFtcFNwZWN0cnVtO1xuICB9LFxuXG4gIFwiemNyXCI6IGZ1bmN0aW9uKGJ1ZmZlclNpemUsIG0pIHtcbiAgICB2YXIgemNyID0gMDtcbiAgICBmb3IgKHZhciBpID0gMDsgaSA8IG0uc2lnbmFsLmxlbmd0aDsgaSsrKSB7XG4gICAgICBpZiAoKG0uc2lnbmFsW2ldID49IDAgJiYgbS5zaWduYWxbaSArIDFdIDwgMCkgfHwgKG0uc2lnbmFsW2ldIDwgMCAmJiBtLnNpZ25hbFtpICsgMV0gPj0gMCkpIHtcbiAgICAgICAgemNyKys7XG4gICAgICB9XG4gICAgfVxuICAgIHJldHVybiB6Y3I7XG4gIH0sXG4gIFxuICBcInBvd2VyU3BlY3RydW1cIjogZnVuY3Rpb24oYnVmZmVyU2l6ZSwgbSkge1xuICAgIHZhciBwb3dlclNwZWN0cnVtID0gbmV3IEZsb2F0MzJBcnJheShtLmFtcFNwZWN0cnVtLmxlbmd0aCk7XG4gICAgZm9yICh2YXIgaSA9IDA7IGkgPCBwb3dlclNwZWN0cnVtLmxlbmd0aDsgaSsrKSB7XG4gICAgICBwb3dlclNwZWN0cnVtW2ldID0gTWF0aC5wb3cobS5hbXBTcGVjdHJ1bVtpXSwgMik7XG4gICAgfVxuICAgIHJldHVybiBwb3dlclNwZWN0cnVtO1xuICB9LFxuICBcbiAgXCJsb3VkbmVzc1wiOiBmdW5jdGlvbihidWZmZXJTaXplLCBtKSB7XG5cbiAgICB2YXIgTlVNX0JBUktfQkFORFMgPSAyNDtcbiAgICB2YXIgc3BlY2lmaWMgPSBuZXcgRmxvYXQzMkFycmF5KE5VTV9CQVJLX0JBTkRTKTtcbiAgICB2YXIgdG90ID0gMDtcbiAgICB2YXIgbm9ybWFsaXNlZFNwZWN0cnVtID0gbS5hbXBTcGVjdHJ1bTtcbiAgICB2YXIgYmJMaW1pdHMgPSBuZXcgSW50MzJBcnJheShOVU1fQkFSS19CQU5EUyArIDEpO1xuXG4gICAgYmJMaW1pdHNbMF0gPSAwO1xuICAgIHZhciBjdXJyZW50QmFuZEVuZCA9IHNlbGYuYmFya1NjYWxlW20uYW1wU3BlY3RydW0ubGVuZ3RoIC0gMV0gLyBOVU1fQkFSS19CQU5EUztcbiAgICB2YXIgY3VycmVudEJhbmQgPSAxO1xuICAgIGZvciAodmFyIGkgPSAwOyBpIDwgbS5hbXBTcGVjdHJ1bS5sZW5ndGg7IGkrKykge1xuICAgICAgd2hpbGUgKHNlbGYuYmFya1NjYWxlW2ldID4gY3VycmVudEJhbmRFbmQpIHtcbiAgICAgICAgYmJMaW1pdHNbY3VycmVudEJhbmQrK10gPSBpO1xuICAgICAgICBjdXJyZW50QmFuZEVuZCA9IGN1cnJlbnRCYW5kICogc2VsZi5iYXJrU2NhbGVbbS5hbXBTcGVjdHJ1bS5sZW5ndGggLSAxXSAvIE5VTV9CQVJLX0JBTkRTO1xuICAgICAgfVxuICAgIH1cblxuICAgIGJiTGltaXRzW05VTV9CQVJLX0JBTkRTXSA9IG0uYW1wU3BlY3RydW0ubGVuZ3RoIC0gMTtcblxuICAgIC8vcHJvY2Vzc1xuXG4gICAgZm9yICh2YXIgaSA9IDA7IGkgPCBOVU1fQkFSS19CQU5EUzsgaSsrKSB7XG4gICAgICB2YXIgc3VtID0gMDtcbiAgICAgIGZvciAodmFyIGogPSBiYkxpbWl0c1tpXTsgaiA8IGJiTGltaXRzW2kgKyAxXTsgaisrKSB7XG5cbiAgICAgICAgc3VtICs9IG5vcm1hbGlzZWRTcGVjdHJ1bVtqXTtcbiAgICAgIH1cbiAgICAgIHNwZWNpZmljW2ldID0gTWF0aC5wb3coc3VtLCAwLjIzKTtcbiAgICB9XG5cbiAgICAvL2dldCB0b3RhbCBsb3VkbmVzc1xuICAgIGZvciAodmFyIGkgPSAwOyBpIDwgc3BlY2lmaWMubGVuZ3RoOyBpKyspIHtcbiAgICAgIHRvdCArPSBzcGVjaWZpY1tpXTtcbiAgICB9XG4gICAgcmV0dXJuIHtcbiAgICAgIFwic3BlY2lmaWNcIjogc3BlY2lmaWMsXG4gICAgICBcInRvdGFsXCI6IHRvdFxuICAgIH07XG4gIH0sXG4gIFxuICBcInBlcmNlcHR1YWxTcHJlYWRcIjogZnVuY3Rpb24oYnVmZmVyU2l6ZSwgbSkge1xuICAgIHZhciBsb3VkbmVzcyA9IG0uZmVhdHVyZUV4dHJhY3RvcnMubG91ZG5lc3MoYnVmZmVyU2l6ZSwgbSk7XG5cbiAgICB2YXIgbWF4ID0gMDtcbiAgICBmb3IgKHZhciBpID0gMDsgaSA8IGxvdWRuZXNzLnNwZWNpZmljLmxlbmd0aDsgaSsrKSB7XG4gICAgICBpZiAobG91ZG5lc3Muc3BlY2lmaWNbaV0gPiBtYXgpIHtcbiAgICAgICAgbWF4ID0gbG91ZG5lc3Muc3BlY2lmaWNbaV07XG4gICAgICB9XG4gICAgfVxuXG4gICAgdmFyIHNwcmVhZCA9IE1hdGgucG93KChsb3VkbmVzcy50b3RhbCAtIG1heCkgLyBsb3VkbmVzcy50b3RhbCwgMik7XG5cbiAgICByZXR1cm4gc3ByZWFkO1xuICB9LFxuICBcbiAgXCJwZXJjZXB0dWFsU2hhcnBuZXNzXCI6IGZ1bmN0aW9uKGJ1ZmZlclNpemUsIG0pIHtcbiAgICB2YXIgbG91ZG5lc3MgPSBtLmZlYXR1cmVFeHRyYWN0b3JzLmxvdWRuZXNzKGJ1ZmZlclNpemUsIG0pO1xuICAgIHZhciBzcGVjID0gbG91ZG5lc3Muc3BlY2lmaWM7XG4gICAgdmFyIG91dHB1dCA9IDA7XG5cbiAgICBmb3IgKHZhciBpID0gMDsgaSA8IHNwZWMubGVuZ3RoOyBpKyspIHtcbiAgICAgIGlmIChpIDwgMTUpIHtcbiAgICAgICAgb3V0cHV0ICs9IChpICsgMSkgKiBzcGVjW2kgKyAxXTtcbiAgICAgIH0gZWxzZSB7XG4gICAgICAgIG91dHB1dCArPSAwLjA2NiAqIE1hdGguZXhwKDAuMTcxICogKGkgKyAxKSk7XG4gICAgICB9XG4gICAgfVxuICAgIG91dHB1dCAqPSAwLjExIC8gbG91ZG5lc3MudG90YWw7XG5cbiAgICByZXR1cm4gb3V0cHV0O1xuICB9LFxuICBcbiAgXCJtZmNjXCI6IGZ1bmN0aW9uKGJ1ZmZlclNpemUsIG0pIHtcbiAgICAvL3VzZWQgdHV0b3JpYWwgZnJvbSBodHRwOi8vcHJhY3RpY2FsY3J5cHRvZ3JhcGh5LmNvbS9taXNjZWxsYW5lb3VzL21hY2hpbmUtbGVhcm5pbmcvZ3VpZGUtbWVsLWZyZXF1ZW5jeS1jZXBzdHJhbC1jb2VmZmljaWVudHMtbWZjY3MvXG4gICAgdmFyIHBvd1NwZWMgPSBtLmZlYXR1cmVFeHRyYWN0b3JzLnBvd2VyU3BlY3RydW0oYnVmZmVyU2l6ZSwgbSk7XG4gICAgdmFyIGZyZXFUb01lbCA9IGZ1bmN0aW9uKGZyZXFWYWx1ZSkge1xuICAgICAgdmFyIG1lbFZhbHVlID0gMTEyNSAqIE1hdGgubG9nKDEgKyAoZnJlcVZhbHVlIC8gNzAwKSk7XG4gICAgICByZXR1cm4gbWVsVmFsdWU7XG4gICAgfTtcbiAgICB2YXIgbWVsVG9GcmVxID0gZnVuY3Rpb24obWVsVmFsdWUpIHtcbiAgICAgIHZhciBmcmVxVmFsdWUgPSA3MDAgKiAoTWF0aC5leHAobWVsVmFsdWUgLyAxMTI1KSAtIDEpO1xuICAgICAgcmV0dXJuIGZyZXFWYWx1ZTtcbiAgICB9O1xuICAgIHZhciBudW1GaWx0ZXJzID0gMjY7IC8vMjYgZmlsdGVycyBpcyBzdGFuZGFyZFxuICAgIHZhciBtZWxWYWx1ZXMgPSBuZXcgRmxvYXQzMkFycmF5KG51bUZpbHRlcnMgKyAyKTsgLy90aGUgKzIgaXMgdGhlIHVwcGVyIGFuZCBsb3dlciBsaW1pdHNcbiAgICB2YXIgbWVsVmFsdWVzSW5GcmVxID0gbmV3IEZsb2F0MzJBcnJheShudW1GaWx0ZXJzICsgMik7XG4gICAgLy9HZW5lcmF0ZSBsaW1pdHMgaW4gSHogLSBmcm9tIDAgdG8gdGhlIG55cXVpc3QuXG4gICAgdmFyIGxvd2VyTGltaXRGcmVxID0gMDtcbiAgICB2YXIgdXBwZXJMaW1pdEZyZXEgPSBhdWRpb0NvbnRleHQuc2FtcGxlUmF0ZSAvIDI7XG4gICAgLy9Db252ZXJ0IHRoZSBsaW1pdHMgdG8gTWVsXG4gICAgdmFyIGxvd2VyTGltaXRNZWwgPSBmcmVxVG9NZWwobG93ZXJMaW1pdEZyZXEpO1xuICAgIHZhciB1cHBlckxpbWl0TWVsID0gZnJlcVRvTWVsKHVwcGVyTGltaXRGcmVxKTtcbiAgICAvL0ZpbmQgdGhlIHJhbmdlXG4gICAgdmFyIHJhbmdlID0gdXBwZXJMaW1pdE1lbCAtIGxvd2VyTGltaXRNZWw7XG4gICAgLy9GaW5kIHRoZSByYW5nZSBhcyBwYXJ0IG9mIHRoZSBsaW5lYXIgaW50ZXJwb2xhdGlvblxuICAgIHZhciB2YWx1ZVRvQWRkID0gcmFuZ2UgLyAobnVtRmlsdGVycyArIDEpO1xuXG4gICAgdmFyIGZmdEJpbnNPZkZyZXEgPSBBcnJheShudW1GaWx0ZXJzICsgMik7XG5cbiAgICBmb3IgKHZhciBpID0gMDsgaSA8IG1lbFZhbHVlcy5sZW5ndGg7IGkrKykge1xuICAgICAgLy9Jbml0aWFsaXNpbmcgdGhlIG1lbCBmcmVxdWVuY2llcyAtIHRoZXkgYXJlIGp1c3QgYSBsaW5lYXIgaW50ZXJwb2xhdGlvbiBiZXR3ZWVuIHRoZSBsb3dlciBhbmQgdXBwZXIgbGltaXRzLlxuICAgICAgbWVsVmFsdWVzW2ldID0gaSAqIHZhbHVlVG9BZGQ7XG4gICAgICAvL0NvbnZlcnQgYmFjayB0byBIelxuICAgICAgbWVsVmFsdWVzSW5GcmVxW2ldID0gbWVsVG9GcmVxKG1lbFZhbHVlc1tpXSk7XG4gICAgICAvL0ZpbmQgdGhlIGNvcnJlc3BvbmRpbmcgYmluc1xuICAgICAgZmZ0Qmluc09mRnJlcVtpXSA9IE1hdGguZmxvb3IoKGJ1ZmZlclNpemUgKyAxKSAqIG1lbFZhbHVlc0luRnJlcVtpXSAvIGF1ZGlvQ29udGV4dC5zYW1wbGVSYXRlKTtcbiAgICB9XG5cbiAgICB2YXIgZmlsdGVyQmFuayA9IEFycmF5KG51bUZpbHRlcnMpO1xuICAgIGZvciAodmFyIGogPSAwOyBqIDwgZmlsdGVyQmFuay5sZW5ndGg7IGorKykge1xuICAgICAgLy9jcmVhdGluZyBhIHR3byBkaW1lbnNpb25hbCBhcnJheSBvZiBzaXplIG51bUZpbHRlcyAqIChidWZmZXJzaXplLzIpKzEgYW5kIHByZS1wb3B1bGF0aW5nIHRoZSBhcnJheXMgd2l0aCAwcy5cbiAgICAgIGZpbHRlckJhbmtbal0gPSBBcnJheS5hcHBseShudWxsLCBuZXcgQXJyYXkoKGJ1ZmZlclNpemUgLyAyKSArIDEpKS5tYXAoTnVtYmVyLnByb3RvdHlwZS52YWx1ZU9mLCAwKTtcbiAgICAgIC8vY3JlYXRpbmcgdGhlIGxvd2VyIGFuZCB1cHBlciBzbG9wZXMgZm9yIGVhY2ggYmluXG4gICAgICBmb3IgKHZhciBpID0gZmZ0Qmluc09mRnJlcVtqXTsgaSA8IGZmdEJpbnNPZkZyZXFbaiArIDFdOyBpKyspIHtcbiAgICAgICAgZmlsdGVyQmFua1tqXVtpXSA9IChpIC0gZmZ0Qmluc09mRnJlcVtqXSkgLyAoZmZ0Qmluc09mRnJlcVtqICsgMV0gLSBmZnRCaW5zT2ZGcmVxW2pdKTtcbiAgICAgIH1cbiAgICAgIGZvciAodmFyIGkgPSBmZnRCaW5zT2ZGcmVxW2ogKyAxXTsgaSA8IGZmdEJpbnNPZkZyZXFbaiArIDJdOyBpKyspIHtcbiAgICAgICAgZmlsdGVyQmFua1tqXVtpXSA9IChmZnRCaW5zT2ZGcmVxW2ogKyAyXSAtIGkpIC8gKGZmdEJpbnNPZkZyZXFbaiArIDJdIC0gZmZ0Qmluc09mRnJlcVtqICsgMV0pO1xuICAgICAgfVxuICAgIH1cblxuICAgIHZhciBsb2dnZWRNZWxCYW5kcyA9IG5ldyBGbG9hdDMyQXJyYXkobnVtRmlsdGVycyk7XG4gICAgZm9yICh2YXIgaSA9IDA7IGkgPCBsb2dnZWRNZWxCYW5kcy5sZW5ndGg7IGkrKykge1xuICAgICAgbG9nZ2VkTWVsQmFuZHNbaV0gPSAwO1xuICAgICAgZm9yICh2YXIgaiA9IDA7IGogPCAoYnVmZmVyU2l6ZSAvIDIpOyBqKyspIHtcbiAgICAgICAgLy9wb2ludCBtdWx0aXBsaWNhdGlvbiBiZXR3ZWVuIHBvd2VyIHNwZWN0cnVtIGFuZCBmaWx0ZXJiYW5rcy5cbiAgICAgICAgZmlsdGVyQmFua1tpXVtqXSA9IGZpbHRlckJhbmtbaV1bal0gKiBwb3dTcGVjW2pdO1xuXG4gICAgICAgIC8vc3VtbWluZyB1cCBhbGwgb2YgdGhlIGNvZWZmaWNpZW50cyBpbnRvIG9uZSBhcnJheVxuICAgICAgICBsb2dnZWRNZWxCYW5kc1tpXSArPSBmaWx0ZXJCYW5rW2ldW2pdO1xuICAgICAgfVxuICAgICAgLy9sb2cgZWFjaCBjb2VmZmljaWVudFxuICAgICAgbG9nZ2VkTWVsQmFuZHNbaV0gPSBNYXRoLmxvZyhsb2dnZWRNZWxCYW5kc1tpXSk7XG4gICAgfVxuXG4gICAgLy9kY3RcbiAgICB2YXIgayA9IE1hdGguUEkgLyBudW1GaWx0ZXJzO1xuICAgIHZhciB3MSA9IDEuMCAvIE1hdGguc3FydChudW1GaWx0ZXJzKTtcbiAgICB2YXIgdzIgPSBNYXRoLnNxcnQoMi4wIC8gbnVtRmlsdGVycyk7XG4gICAgdmFyIG51bUNvZWZmcyA9IDEzO1xuICAgIHZhciBkY3RNYXRyaXggPSBuZXcgRmxvYXQzMkFycmF5KG51bUNvZWZmcyAqIG51bUZpbHRlcnMpO1xuXG4gICAgZm9yICh2YXIgaSA9IDA7IGkgPCBudW1Db2VmZnM7IGkrKykge1xuICAgICAgZm9yICh2YXIgaiA9IDA7IGogPCBudW1GaWx0ZXJzOyBqKyspIHtcbiAgICAgICAgdmFyIGlkeCA9IGkgKyAoaiAqIG51bUNvZWZmcyk7XG4gICAgICAgIGlmIChpID09IDApIHtcbiAgICAgICAgICBkY3RNYXRyaXhbaWR4XSA9IHcxICogTWF0aC5jb3MoayAqIChpICsgMSkgKiAoaiArIDAuNSkpO1xuICAgICAgICB9IGVsc2Uge1xuICAgICAgICAgIGRjdE1hdHJpeFtpZHhdID0gdzIgKiBNYXRoLmNvcyhrICogKGkgKyAxKSAqIChqICsgMC41KSk7XG4gICAgICAgIH1cbiAgICAgIH1cbiAgICB9XG5cbiAgICB2YXIgbWZjY3MgPSBuZXcgRmxvYXQzMkFycmF5KG51bUNvZWZmcyk7XG4gICAgZm9yICh2YXIgayA9IDA7IGsgPCBudW1Db2VmZnM7IGsrKykge1xuICAgICAgdmFyIHYgPSAwO1xuICAgICAgZm9yICh2YXIgbiA9IDA7IG4gPCBudW1GaWx0ZXJzOyBuKyspIHtcbiAgICAgICAgdmFyIGlkeCA9IGsgKyAobiAqIG51bUNvZWZmcyk7XG4gICAgICAgIHYgKz0gKGRjdE1hdHJpeFtpZHhdICogbG9nZ2VkTWVsQmFuZHNbbl0pO1xuICAgICAgfVxuICAgICAgbWZjY3Nba10gPSB2IC8gbnVtQ29lZmZzO1xuICAgIH1cbiAgICByZXR1cm4gbWZjY3M7XG4gIH1cbn07IiwiLy8gTWV5ZGEgSmF2YXNjcmlwdCBEU1AgbGlicmFyeVxuXG52YXIgQ29tcGxleEFycmF5ID0gcmVxdWlyZSgnLi4vbGliL2pzZmZ0L2NvbXBsZXhfYXJyYXknKS5Db21wbGV4QXJyYXk7XG4vLyBtb2RpZmllcyBDb21wbGV4QXJyYXlcbnZhciBmZnQgPSByZXF1aXJlKCcuLi9saWIvanNmZnQvZmZ0Jyk7XG5cbm1vZHVsZS5leHBvcnRzID0gZnVuY3Rpb24oYXVkaW9Db250ZXh0LHNyYyxidWZTaXplLGNhbGxiYWNrKXtcblx0XG5cdC8vSSBhbSBteXNlbGZcblx0dmFyIHNlbGYgPSB0aGlzO1xuXHRcblx0c2VsZi5mZWF0dXJlRXh0cmFjdG9ycyA9IHJlcXVpcmUoJy4vZXh0cmFjdG9ycycpO1xuXG5cdC8vZGVmYXVsdCBidWZmZXIgc2l6ZVxuXHR2YXIgYnVmZmVyU2l6ZSA9IGJ1ZlNpemUgPyBidWZTaXplIDogMjU2O1xuXG5cdC8vaW5pdGlhbCBzb3VyY2Vcblx0dmFyIHNvdXJjZSA9IHNyYztcblxuXHQvL2NhbGxiYWNrIGNvbnRyb2xsZXJzXG5cdHZhciBFWFRSQUNUSU9OX1NUQVJURUQgPSBmYWxzZTtcblx0dmFyIF9mZWF0dXJlc1RvRXh0cmFjdDtcblxuXHQvL3V0aWxpdGllc1xuXHR2YXIgwrUgPSBmdW5jdGlvbihpLCBhbXBsaXR1ZGVTcGVjdCl7XG5cdFx0dmFyIG51bWVyYXRvciA9IDA7XG5cdFx0dmFyIGRlbm9taW5hdG9yID0gMDtcblx0XHRmb3IodmFyIGsgPSAwOyBrIDwgYW1wbGl0dWRlU3BlY3QubGVuZ3RoOyBrKyspe1xuXHRcdFx0bnVtZXJhdG9yICs9IE1hdGgucG93KGssaSkqTWF0aC5hYnMoYW1wbGl0dWRlU3BlY3Rba10pO1xuXHRcdFx0ZGVub21pbmF0b3IgKz0gYW1wbGl0dWRlU3BlY3Rba107XG5cdFx0fVxuXHRcdHJldHVybiBudW1lcmF0b3IvZGVub21pbmF0b3I7XG5cdH07XG5cblx0dmFyIGlzUG93ZXJPZlR3byA9IGZ1bmN0aW9uKG51bSkge1xuXHRcdHdoaWxlICgoKG51bSAlIDIpID09IDApICYmIG51bSA+IDEpIHtcblx0XHRcdG51bSAvPSAyO1xuXHRcdH1cblx0XHRyZXR1cm4gKG51bSA9PSAxKTtcblx0fTtcblxuXHQvL2luaXRpbGl6ZSBiYXJrIHNjYWxlICh0byBiZSB1c2VkIGluIG1vc3QgcGVyY2VwdHVhbCBmZWF0dXJlcykuXG5cdHNlbGYuYmFya1NjYWxlID0gbmV3IEZsb2F0MzJBcnJheShidWZTaXplKTtcblxuXHRmb3IodmFyIGkgPSAwOyBpIDwgc2VsZi5iYXJrU2NhbGUubGVuZ3RoOyBpKyspe1xuXHRcdHNlbGYuYmFya1NjYWxlW2ldID0gaSphdWRpb0NvbnRleHQuc2FtcGxlUmF0ZS8oYnVmU2l6ZSk7XG5cdFx0c2VsZi5iYXJrU2NhbGVbaV0gPSAxMypNYXRoLmF0YW4oc2VsZi5iYXJrU2NhbGVbaV0vMTMxNS44KSArIDMuNSogTWF0aC5hdGFuKE1hdGgucG93KChzZWxmLmJhcmtTY2FsZVtpXS83NTE4KSwyKSk7XG5cdH1cblxuXHQvL1dJTkRPV0lOR1xuXHQvL3NldCBkZWZhdWx0XG5cdHNlbGYud2luZG93aW5nRnVuY3Rpb24gPSBcImhhbm5pbmdcIjtcblxuXHQvL2NyZWF0ZSB3aW5kb3dzXG5cdHNlbGYuaGFubmluZyA9IG5ldyBGbG9hdDMyQXJyYXkoYnVmU2l6ZSk7XG5cdGZvciAodmFyIGkgPSAwOyBpIDwgYnVmU2l6ZTsgaSsrKSB7XG5cdFx0Ly9BY2NvcmRpbmcgdG8gdGhlIFIgZG9jdW1lbnRhdGlvbiBodHRwOi8vcmdtLm9nYWxhYi5uZXQvUkdNL1JfcmRmaWxlP2Y9R0VORUFyZWFkL21hbi9oYW5uaW5nLndpbmRvdy5SZCZkPVJfQ0Ncblx0XHRzZWxmLmhhbm5pbmdbaV0gPSAwLjUgLSAwLjUqTWF0aC5jb3MoMipNYXRoLlBJKmkvKGJ1ZlNpemUtMSkpO1xuXHR9XG5cblx0c2VsZi5oYW1taW5nID0gbmV3IEZsb2F0MzJBcnJheShidWZTaXplKTtcblx0Zm9yICh2YXIgaSA9IDA7IGkgPCBidWZTaXplOyBpKyspIHtcblx0XHQvL0FjY29yZGluZyB0byBodHRwOi8vdWsubWF0aHdvcmtzLmNvbS9oZWxwL3NpZ25hbC9yZWYvaGFtbWluZy5odG1sXG5cdFx0c2VsZi5oYW1taW5nW2ldID0gMC41NCAtIDAuNDYqTWF0aC5jb3MoMipNYXRoLlBJKihpL2J1ZlNpemUtMSkpO1xuXHR9XG5cblx0Ly9VTkZJTklTSEVEIC0gYmxhY2ttYW4gd2luZG93IGltcGxlbWVudGF0aW9uXG5cblx0LypzZWxmLmJsYWNrbWFuID0gbmV3IEZsb2F0MzJBcnJheShidWZTaXplKTtcblx0Ly9BY2NvcmRpbmcgdG8gaHR0cDovL3VrLm1hdGh3b3Jrcy5jb20vaGVscC9zaWduYWwvcmVmL2JsYWNrbWFuLmh0bWxcblx0Ly9maXJzdCBoYWxmIG9mIHRoZSB3aW5kb3dcblx0Zm9yICh2YXIgaSA9IDA7IGkgPCAoYnVmU2l6ZSAlIDIpID8gKGJ1ZlNpemUrMSkvMiA6IGJ1ZlNpemUvMjsgaSsrKSB7XG5cdFx0c2VsZi5ibGFja21hbltpXSA9IDAuNDIgLSAwLjUqTWF0aC5jb3MoMipNYXRoLlBJKmkvKGJ1ZlNpemUtMSkpICsgMC4wOCpNYXRoLmNvcyg0Kk1hdGguUEkqaS8oYnVmU2l6ZS0xKSk7XG5cdH1cblx0Ly9zZWNvbmQgaGFsZiBvZiB0aGUgd2luZG93XG5cdGZvciAodmFyIGkgPSBidWZTaXplLzI7IGkgPiAwOyBpLS0pIHtcblx0XHRzZWxmLmJsYWNrbWFuW2J1ZlNpemUgLSBpXSA9IHNlbGYuYmxhY2ttYW5baV07XG5cdH0qL1xuXG5cdHNlbGYud2luZG93aW5nID0gZnVuY3Rpb24oc2lnLCB0eXBlKXtcblx0XHR2YXIgd2luZG93ZWQgPSBuZXcgRmxvYXQzMkFycmF5KHNpZy5sZW5ndGgpO1xuXHRcdHZhciBpICxsZW4gPSBzaWcubGVuZ3RoO1xuXG5cdFx0aWYgKHR5cGUgPT0gXCJoYW5uaW5nXCIpIHtcblx0XHRcdGZvciAoaSA9IDA7IGkgPCBsZW47IGkrKykge1xuXHRcdFx0XHR3aW5kb3dlZFtpXSA9IHNpZ1tpXSpzZWxmLmhhbm5pbmdbaV07XG5cdFx0XHR9XG5cdFx0fVxuXHRcdGVsc2UgaWYgKHR5cGUgPT0gXCJoYW1taW5nXCIpIHtcblx0XHRcdGZvciAoaSA9IDA7IGkgPCBsZW47IGkrKykge1xuXHRcdFx0XHR3aW5kb3dlZFtpXSA9IHNpZ1tpXSpzZWxmLmhhbW1pbmdbaV07XG5cdFx0XHR9XG5cdFx0fVxuXHRcdGVsc2UgaWYgKHR5cGUgPT0gXCJibGFja21hblwiKSB7XG5cdFx0XHRmb3IgKGkgPSAwOyBpIDwgbGVuOyBpKyspIHtcblx0XHRcdFx0d2luZG93ZWRbaV0gPSBzaWdbaV0qc2VsZi5ibGFja21hbltpXTtcblx0XHRcdH1cblx0XHR9XG5cblx0XHRyZXR1cm4gd2luZG93ZWQ7XG5cdH07XG5cblx0Ly9zb3VyY2Ugc2V0dGVyIG1ldGhvZFxuXHRzZWxmLnNldFNvdXJjZSA9IGZ1bmN0aW9uKF9zcmMpIHtcblx0XHRzb3VyY2UgPSBfc3JjO1xuXHRcdHNvdXJjZS5jb25uZWN0KHdpbmRvdy5zcG4pO1xuXHR9O1xuXG5cblxuXHRpZiAoaXNQb3dlck9mVHdvKGJ1ZmZlclNpemUpICYmIGF1ZGlvQ29udGV4dCkge1xuXHRcdFx0c2VsZi5mZWF0dXJlSW5mbyA9IHtcblx0XHRcdFx0XCJidWZmZXJcIjoge1xuXHRcdFx0XHRcdFwidHlwZVwiOiBcImFycmF5XCJcblx0XHRcdFx0fSxcblx0XHRcdFx0XCJybXNcIjoge1xuXHRcdFx0XHRcdFwidHlwZVwiOiBcIm51bWJlclwiXG5cdFx0XHRcdH0sXG5cdFx0XHRcdFwiZW5lcmd5XCI6IHtcblx0XHRcdFx0XHRcInR5cGVcIjogXCJudW1iZXJcIlxuXHRcdFx0XHR9LFxuXHRcdFx0XHRcInpjclwiOiB7XG5cdFx0XHRcdFx0XCJ0eXBlXCI6IFwibnVtYmVyXCJcblx0XHRcdFx0fSxcblx0XHRcdFx0XCJjb21wbGV4U3BlY3RydW1cIjoge1xuXHRcdFx0XHRcdFwidHlwZVwiOiBcIm11bHRpcGxlQXJyYXlzXCIsXG5cdFx0XHRcdFx0XCJhcnJheU5hbWVzXCI6IHtcblx0XHRcdFx0XHRcdFwiMVwiOiBcInJlYWxcIixcblx0XHRcdFx0XHRcdFwiMlwiOiBcImltYWdcIlxuXHRcdFx0XHRcdH1cblx0XHRcdFx0fSxcblx0XHRcdFx0XCJhbXBsaXR1ZGVTcGVjdHJ1bVwiOiB7XG5cdFx0XHRcdFx0XCJ0eXBlXCI6IFwiYXJyYXlcIlxuXHRcdFx0XHR9LFxuXHRcdFx0XHRcInBvd2VyU3BlY3RydW1cIjoge1xuXHRcdFx0XHRcdFwidHlwZVwiOiBcImFycmF5XCJcblx0XHRcdFx0fSxcblx0XHRcdFx0XCJzcGVjdHJhbENlbnRyb2lkXCI6IHtcblx0XHRcdFx0XHRcInR5cGVcIjogXCJudW1iZXJcIlxuXHRcdFx0XHR9LFxuXHRcdFx0XHRcInNwZWN0cmFsRmxhdG5lc3NcIjoge1xuXHRcdFx0XHRcdFwidHlwZVwiOiBcIm51bWJlclwiXG5cdFx0XHRcdH0sXG5cdFx0XHRcdFwic3BlY3RyYWxTbG9wZVwiOiB7XG5cdFx0XHRcdFx0XCJ0eXBlXCI6IFwibnVtYmVyXCJcblx0XHRcdFx0fSxcblx0XHRcdFx0XCJzcGVjdHJhbFJvbGxvZmZcIjoge1xuXHRcdFx0XHRcdFwidHlwZVwiOiBcIm51bWJlclwiXG5cdFx0XHRcdH0sXG5cdFx0XHRcdFwic3BlY3RyYWxTcHJlYWRcIjoge1xuXHRcdFx0XHRcdFwidHlwZVwiOiBcIm51bWJlclwiXG5cdFx0XHRcdH0sXG5cdFx0XHRcdFwic3BlY3RyYWxTa2V3bmVzc1wiOiB7XG5cdFx0XHRcdFx0XCJ0eXBlXCI6IFwibnVtYmVyXCJcblx0XHRcdFx0fSxcblx0XHRcdFx0XCJzcGVjdHJhbEt1cnRvc2lzXCI6IHtcblx0XHRcdFx0XHRcInR5cGVcIjogXCJudW1iZXJcIlxuXHRcdFx0XHR9LFxuXHRcdFx0XHRcImxvdWRuZXNzXCI6IHtcblx0XHRcdFx0XHRcInR5cGVcIjogXCJtdWx0aXBsZUFycmF5c1wiLFxuXHRcdFx0XHRcdFwiYXJyYXlOYW1lc1wiOiB7XG5cdFx0XHRcdFx0XHRcIjFcIjogXCJ0b3RhbFwiLFxuXHRcdFx0XHRcdFx0XCIyXCI6IFwic3BlY2lmaWNcIlxuXHRcdFx0XHRcdH1cblx0XHRcdFx0fSxcblx0XHRcdFx0XCJwZXJjZXB0dWFsU3ByZWFkXCI6IHtcblx0XHRcdFx0XHRcInR5cGVcIjogXCJudW1iZXJcIlxuXHRcdFx0XHR9LFxuXHRcdFx0XHRcInBlcmNlcHR1YWxTaGFycG5lc3NcIjoge1xuXHRcdFx0XHRcdFwidHlwZVwiOiBcIm51bWJlclwiXG5cdFx0XHRcdH0sXG5cdFx0XHRcdFwibWZjY1wiOiB7XG5cdFx0XHRcdFx0XCJ0eXBlXCI6IFwiYXJyYXlcIlxuXHRcdFx0XHR9XG5cdFx0XHR9O1xuXG5cdFx0XHQvL2NyZWF0ZSBjb21wbGV4YXJyYXkgdG8gaG9sZCB0aGUgc3BlY3RydW1cblx0XHRcdHZhciBkYXRhID0gbmV3IENvbXBsZXhBcnJheShidWZmZXJTaXplKTtcblxuXHRcdFx0Ly9jcmVhdGUgbm9kZXNcblx0XHRcdHdpbmRvdy5zcG4gPSBhdWRpb0NvbnRleHQuY3JlYXRlU2NyaXB0UHJvY2Vzc29yKGJ1ZmZlclNpemUsMSwxKTtcblx0XHRcdHNwbi5jb25uZWN0KGF1ZGlvQ29udGV4dC5kZXN0aW5hdGlvbik7XG5cblx0XHRcdHdpbmRvdy5zcG4ub25hdWRpb3Byb2Nlc3MgPSBmdW5jdGlvbihlKSB7XG5cdFx0XHRcdC8vdGhpcyBpcyB0byBvYnRhaW4gdGhlIGN1cnJlbnQgYW1wbGl0dWRlIHNwZWN0cnVtXG5cdFx0XHRcdHZhciBpbnB1dERhdGEgPSBlLmlucHV0QnVmZmVyLmdldENoYW5uZWxEYXRhKDApO1xuXHRcdFx0XHRzZWxmLnNpZ25hbCA9IGlucHV0RGF0YTtcblx0XHRcdFx0dmFyIHdpbmRvd2VkU2lnbmFsID0gc2VsZi53aW5kb3dpbmcoc2VsZi5zaWduYWwsIHNlbGYud2luZG93aW5nRnVuY3Rpb24pO1xuXG5cdFx0XHRcdC8vbWFwIHRpbWUgZG9tYWluXG5cdFx0XHRcdGRhdGEubWFwKGZ1bmN0aW9uKHZhbHVlLCBpLCBuKSB7XG5cdFx0XHRcdFx0dmFsdWUucmVhbCA9IHdpbmRvd2VkU2lnbmFsW2ldO1xuXHRcdFx0XHR9KTtcblx0XHRcdFx0Ly90cmFuc2Zvcm1cblx0XHRcdFx0dmFyIHNwZWMgPSBkYXRhLkZGVCgpO1xuXHRcdFx0XHQvL2Fzc2lnbiB0byBtZXlkYVxuXHRcdFx0XHRzZWxmLmNvbXBsZXhTcGVjdHJ1bSA9IHNwZWM7XG5cdFx0XHRcdHNlbGYuYW1wU3BlY3RydW0gPSBuZXcgRmxvYXQzMkFycmF5KGJ1ZmZlclNpemUvMik7XG5cdFx0XHRcdC8vY2FsY3VsYXRlIGFtcGxpdHVkZVxuXHRcdFx0XHRmb3IgKHZhciBpID0gMDsgaSA8IGJ1ZmZlclNpemUvMjsgaSsrKSB7XG5cdFx0XHRcdFx0c2VsZi5hbXBTcGVjdHJ1bVtpXSA9IE1hdGguc3FydChNYXRoLnBvdyhzcGVjLnJlYWxbaV0sMikgKyBNYXRoLnBvdyhzcGVjLmltYWdbaV0sMikpO1xuXG5cdFx0XHRcdH1cblx0XHRcdFx0Ly9jYWxsIGNhbGxiYWNrIGlmIGFwcGxpY2FibGVcblx0XHRcdFx0aWYgKHR5cGVvZiBjYWxsYmFjayA9PT0gXCJmdW5jdGlvblwiICYmIEVYVFJBQ1RJT05fU1RBUlRFRCkge1xuXHRcdFx0XHRcdGNhbGxiYWNrKHNlbGYuZ2V0KF9mZWF0dXJlc1RvRXh0cmFjdCkpO1xuXHRcdFx0XHR9XG5cblx0XHRcdH07XG5cblx0XHRcdHNlbGYuc3RhcnQgPSBmdW5jdGlvbihmZWF0dXJlcykge1xuXHRcdFx0XHRfZmVhdHVyZXNUb0V4dHJhY3QgPSBmZWF0dXJlcztcblx0XHRcdFx0RVhUUkFDVElPTl9TVEFSVEVEID0gdHJ1ZTtcblx0XHRcdH07XG5cblx0XHRcdHNlbGYuc3RvcCA9IGZ1bmN0aW9uKCkge1xuXHRcdFx0XHRFWFRSQUNUSU9OX1NUQVJURUQgPSBmYWxzZTtcblx0XHRcdH07XG5cblx0XHRcdHNlbGYuYXVkaW9Db250ZXh0ID0gYXVkaW9Db250ZXh0O1xuXG5cdFx0XHRzZWxmLmdldCA9IGZ1bmN0aW9uKGZlYXR1cmUpIHtcblx0XHRcdFx0aWYodHlwZW9mIGZlYXR1cmUgPT09IFwib2JqZWN0XCIpe1xuXHRcdFx0XHRcdHZhciByZXN1bHRzID0ge307XG5cdFx0XHRcdFx0Zm9yICh2YXIgeCA9IDA7IHggPCBmZWF0dXJlLmxlbmd0aDsgeCsrKXtcblx0XHRcdFx0XHRcdHRyeXtcblx0XHRcdFx0XHRcdFx0cmVzdWx0c1tmZWF0dXJlW3hdXSA9IChzZWxmLmZlYXR1cmVFeHRyYWN0b3JzW2ZlYXR1cmVbeF1dKGJ1ZmZlclNpemUsIHNlbGYpKTtcblx0XHRcdFx0XHRcdH0gY2F0Y2ggKGUpe1xuXHRcdFx0XHRcdFx0XHRjb25zb2xlLmVycm9yKGUpO1xuXHRcdFx0XHRcdFx0fVxuXHRcdFx0XHRcdH1cblx0XHRcdFx0XHRyZXR1cm4gcmVzdWx0cztcblx0XHRcdFx0fVxuXHRcdFx0XHRlbHNlIGlmICh0eXBlb2YgZmVhdHVyZSA9PT0gXCJzdHJpbmdcIil7XG5cdFx0XHRcdFx0cmV0dXJuIHNlbGYuZmVhdHVyZUV4dHJhY3RvcnNbZmVhdHVyZV0oYnVmZmVyU2l6ZSwgc2VsZik7XG5cdFx0XHRcdH1cblx0XHRcdFx0ZWxzZXtcblx0XHRcdFx0XHR0aHJvdyBcIkludmFsaWQgRmVhdHVyZSBGb3JtYXRcIjtcblx0XHRcdFx0fVxuXHRcdFx0fTtcblx0XHRcdHNvdXJjZS5jb25uZWN0KHdpbmRvdy5zcG4sIDAsIDApO1xuXHRcdFx0cmV0dXJuIHNlbGY7XG5cdH1cblx0ZWxzZSB7XG5cdFx0Ly9oYW5kbGUgZXJyb3JzXG5cdFx0aWYgKHR5cGVvZiBhdWRpb0NvbnRleHQgPT0gXCJ1bmRlZmluZWRcIikge1xuXHRcdFx0dGhyb3cgXCJBdWRpb0NvbnRleHQgd2Fzbid0IHNwZWNpZmllZDogTWV5ZGEgd2lsbCBub3QgcnVuLlwiO1xuXHRcdH1cblx0XHRlbHNlIHtcblx0XHRcdHRocm93IFwiQnVmZmVyIHNpemUgaXMgbm90IGEgcG93ZXIgb2YgdHdvOiBNZXlkYSB3aWxsIG5vdCBydW4uXCI7XG5cdFx0fVxuXHR9XG59OyJdfQ==

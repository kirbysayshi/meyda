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
module.exports = function(bufferSize, m){
  return m.ampSpectrum;
};
},{}],4:[function(require,module,exports){
module.exports = function(bufferSize, m) {
  return m.signal;
};
},{}],5:[function(require,module,exports){
module.exports = function(bufferSize, m) {
  return m.complexSpectrum;
};
},{}],6:[function(require,module,exports){
module.exports = function(bufferSize, m) {
  var energy = 0;
  for (var i = 0; i < m.signal.length; i++) {
    energy += Math.pow(Math.abs(m.signal[i]), 2);
  }
  return energy;
};
},{}],7:[function(require,module,exports){
module.exports = {
  buffer: require('./buffer'),
  rms: require('./rms'),
  energy: require('./energy'),
  complexSpectrum: require('./complexSpectrum'),
  spectralSlope: require('./spectralSlope'),
  spectralCentroid: require('./spectralCentroid'),
  spectralRolloff: require('./spectralRolloff'),
  spectralFlatness: require('./spectralFlatness'),
  spectralSpread: require('./spectralSpread'),
  spectralSkewness: require('./spectralSkewness'),
  spectralKurtosis: require('./spectralKurtosis'),
  amplitudeSpectrum: require('./amplitudeSpectrum'),
  zcr: require('./zcr'),
  powerSpectrum: require('./powerSpectrum'),
  loudness: require('./loudness'),
  perceptualSpread: require('./perceptualSpread'),
  perceptualSharpness: require('./perceptualSharpness'),
  mfcc: require('./mfcc')
};
},{"./amplitudeSpectrum":3,"./buffer":4,"./complexSpectrum":5,"./energy":6,"./loudness":8,"./mfcc":9,"./perceptualSharpness":10,"./perceptualSpread":11,"./powerSpectrum":12,"./rms":13,"./spectralCentroid":14,"./spectralFlatness":15,"./spectralKurtosis":16,"./spectralRolloff":17,"./spectralSkewness":18,"./spectralSlope":19,"./spectralSpread":20,"./zcr":21}],8:[function(require,module,exports){
module.exports = function(bufferSize, m) {
  var barkScale = new Float32Array(m.ampSpectrum.length);
  var NUM_BARK_BANDS = 24;
  var specific = new Float32Array(NUM_BARK_BANDS);
  var tot = 0;
  var normalisedSpectrum = m.ampSpectrum;
  var bbLimits = new Int32Array(NUM_BARK_BANDS + 1);

  for (var i = 0; i < barkScale.length; i++) {
    barkScale[i] = i * m.audioContext.sampleRate / (bufferSize);
    barkScale[i] = 13 * Math.atan(barkScale[i] / 1315.8) + 3.5 * Math.atan(Math.pow((barkScale[i] / 7518), 2));
  }


  bbLimits[0] = 0;
  var currentBandEnd = barkScale[m.ampSpectrum.length - 1] / NUM_BARK_BANDS;
  var currentBand = 1;
  for (var i = 0; i < m.ampSpectrum.length; i++) {
    while (barkScale[i] > currentBandEnd) {
      bbLimits[currentBand++] = i;
      currentBandEnd = currentBand * barkScale[m.ampSpectrum.length - 1] / NUM_BARK_BANDS;
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
};
},{}],9:[function(require,module,exports){
//used tutorial from http://practicalcryptography.com/miscellaneous/machine-learning/guide-mel-frequency-cepstral-coefficients-mfccs/

var powerSpectrum = require('./powerSpectrum');

module.exports = function(bufferSize, m) {
  var powSpec = powerSpectrum(bufferSize, m);
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
};
},{"./powerSpectrum":12}],10:[function(require,module,exports){
module.exports = function(bufferSize, m) {
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
};
},{}],11:[function(require,module,exports){
module.exports = function(bufferSize, m) {
  var loudness = m.featureExtractors.loudness(bufferSize, m);

  var max = 0;
  for (var i = 0; i < loudness.specific.length; i++) {
    if (loudness.specific[i] > max) {
      max = loudness.specific[i];
    }
  }

  var spread = Math.pow((loudness.total - max) / loudness.total, 2);

  return spread;
};
},{}],12:[function(require,module,exports){
module.exports = function(bufferSize, m) {
  var powerSpectrum = new Float32Array(m.ampSpectrum.length);
  for (var i = 0; i < powerSpectrum.length; i++) {
    powerSpectrum[i] = Math.pow(m.ampSpectrum[i], 2);
  }
  return powerSpectrum;
};
},{}],13:[function(require,module,exports){
module.exports = function(bufferSize, m) {

  var rms = 0;
  for (var i = 0; i < m.signal.length; i++) {
    rms += Math.pow(m.signal[i], 2);
  }
  rms = rms / m.signal.length;
  rms = Math.sqrt(rms);

  return rms;
};

},{}],14:[function(require,module,exports){
module.exports = function(bufferSize, m) {
  return µ(1, m.ampSpectrum);
};
},{}],15:[function(require,module,exports){
module.exports = function(bufferSize, m) {
  var ampspec = m.ampSpectrum;
  var numerator = 0;
  var denominator = 0;
  for (var i = 0; i < ampspec.length; i++) {
    numerator += Math.log(ampspec[i]);
    denominator += ampspec[i];
  }
  return Math.exp(numerator / ampspec.length) * ampspec.length / denominator;
};
},{}],16:[function(require,module,exports){
module.exports = function(bufferSize, m) {
  var ampspec = m.ampSpectrum;
  var µ1 = µ(1, ampspec);
  var µ2 = µ(2, ampspec);
  var µ3 = µ(3, ampspec);
  var µ4 = µ(4, ampspec);
  var numerator = -3 * Math.pow(µ1, 4) + 6 * µ1 * µ2 - 4 * µ1 * µ3 + µ4;
  var denominator = Math.pow(Math.sqrt(µ2 - Math.pow(µ1, 2)), 4);
  return numerator / denominator;
};
},{}],17:[function(require,module,exports){
module.exports = function(bufferSize, m) {
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
};
},{}],18:[function(require,module,exports){
module.exports = function(bufferSize, m, spectrum) {
  var ampspec = m.ampSpectrum;
  var µ1 = µ(1, ampspec);
  var µ2 = µ(2, ampspec);
  var µ3 = µ(3, ampspec);
  var numerator = 2 * Math.pow(µ1, 3) - 3 * µ1 * µ2 + µ3;
  var denominator = Math.pow(Math.sqrt(µ2 - Math.pow(µ1, 2)), 3);
  return numerator / denominator;
};
},{}],19:[function(require,module,exports){
module.exports = function(bufferSize, m) {
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
};
},{}],20:[function(require,module,exports){
module.exports = function(bufferSize, m) {
  var ampspec = m.ampSpectrum;
  return Math.sqrt(µ(2, ampspec) - Math.pow(µ(1, ampspec), 2));
};
},{}],21:[function(require,module,exports){
module.exports = function(bufferSize, m) {
  var zcr = 0;
  for (var i = 0; i < m.signal.length; i++) {
    if ((m.signal[i] >= 0 && m.signal[i + 1] < 0) || (m.signal[i] < 0 && m.signal[i + 1] >= 0)) {
      zcr++;
    }
  }
  return zcr;
};
},{}],22:[function(require,module,exports){
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
			
			//transform
			var spec = data.FFT();
			//assign to meyda
			self.complexSpectrum = spec;
			self.ampSpectrum = new Float32Array(bufferSize/2);

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
},{"../lib/jsfft/complex_array":1,"../lib/jsfft/fft":2,"./extractors":7}]},{},[22])(22)
});
//# sourceMappingURL=data:application/json;base64,eyJ2ZXJzaW9uIjozLCJzb3VyY2VzIjpbIi4uLy4uLy4uLy4uL3Vzci9sb2NhbC9saWIvbm9kZV9tb2R1bGVzL2Jyb3dzZXJpZnkvbm9kZV9tb2R1bGVzL2Jyb3dzZXItcGFjay9fcHJlbHVkZS5qcyIsImxpYi9qc2ZmdC9jb21wbGV4X2FycmF5LmpzIiwibGliL2pzZmZ0L2ZmdC5qcyIsInNyYy9leHRyYWN0b3JzL2FtcGxpdHVkZVNwZWN0cnVtLmpzIiwic3JjL2V4dHJhY3RvcnMvYnVmZmVyLmpzIiwic3JjL2V4dHJhY3RvcnMvY29tcGxleFNwZWN0cnVtLmpzIiwic3JjL2V4dHJhY3RvcnMvZW5lcmd5LmpzIiwic3JjL2V4dHJhY3RvcnMvaW5kZXguanMiLCJzcmMvZXh0cmFjdG9ycy9sb3VkbmVzcy5qcyIsInNyYy9leHRyYWN0b3JzL21mY2MuanMiLCJzcmMvZXh0cmFjdG9ycy9wZXJjZXB0dWFsU2hhcnBuZXNzLmpzIiwic3JjL2V4dHJhY3RvcnMvcGVyY2VwdHVhbFNwcmVhZC5qcyIsInNyYy9leHRyYWN0b3JzL3Bvd2VyU3BlY3RydW0uanMiLCJzcmMvZXh0cmFjdG9ycy9ybXMuanMiLCJzcmMvZXh0cmFjdG9ycy9zcGVjdHJhbENlbnRyb2lkLmpzIiwic3JjL2V4dHJhY3RvcnMvc3BlY3RyYWxGbGF0bmVzcy5qcyIsInNyYy9leHRyYWN0b3JzL3NwZWN0cmFsS3VydG9zaXMuanMiLCJzcmMvZXh0cmFjdG9ycy9zcGVjdHJhbFJvbGxvZmYuanMiLCJzcmMvZXh0cmFjdG9ycy9zcGVjdHJhbFNrZXduZXNzLmpzIiwic3JjL2V4dHJhY3RvcnMvc3BlY3RyYWxTbG9wZS5qcyIsInNyYy9leHRyYWN0b3JzL3NwZWN0cmFsU3ByZWFkLmpzIiwic3JjL2V4dHJhY3RvcnMvemNyLmpzIiwic3JjL21leWRhLmpzIl0sIm5hbWVzIjpbXSwibWFwcGluZ3MiOiJBQUFBO0FDQUE7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBOztBQ3BIQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTs7QUNqT0E7QUFDQTtBQUNBOztBQ0ZBO0FBQ0E7QUFDQTs7QUNGQTtBQUNBO0FBQ0E7O0FDRkE7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7O0FDTkE7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTs7QUNuQkE7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7O0FDN0NBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTs7QUMvRkE7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7O0FDZkE7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTs7QUNiQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTs7QUNOQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7O0FDWEE7QUFDQTtBQUNBOztBQ0ZBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBOztBQ1RBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBOztBQ1RBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBOztBQ2ZBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTs7QUNSQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7O0FDakJBO0FBQ0E7QUFDQTtBQUNBOztBQ0hBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTs7QUNSQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQSIsImZpbGUiOiJnZW5lcmF0ZWQuanMiLCJzb3VyY2VSb290IjoiIiwic291cmNlc0NvbnRlbnQiOlsiKGZ1bmN0aW9uIGUodCxuLHIpe2Z1bmN0aW9uIHMobyx1KXtpZighbltvXSl7aWYoIXRbb10pe3ZhciBhPXR5cGVvZiByZXF1aXJlPT1cImZ1bmN0aW9uXCImJnJlcXVpcmU7aWYoIXUmJmEpcmV0dXJuIGEobywhMCk7aWYoaSlyZXR1cm4gaShvLCEwKTt2YXIgZj1uZXcgRXJyb3IoXCJDYW5ub3QgZmluZCBtb2R1bGUgJ1wiK28rXCInXCIpO3Rocm93IGYuY29kZT1cIk1PRFVMRV9OT1RfRk9VTkRcIixmfXZhciBsPW5bb109e2V4cG9ydHM6e319O3Rbb11bMF0uY2FsbChsLmV4cG9ydHMsZnVuY3Rpb24oZSl7dmFyIG49dFtvXVsxXVtlXTtyZXR1cm4gcyhuP246ZSl9LGwsbC5leHBvcnRzLGUsdCxuLHIpfXJldHVybiBuW29dLmV4cG9ydHN9dmFyIGk9dHlwZW9mIHJlcXVpcmU9PVwiZnVuY3Rpb25cIiYmcmVxdWlyZTtmb3IodmFyIG89MDtvPHIubGVuZ3RoO28rKylzKHJbb10pO3JldHVybiBzfSkiLCIndXNlIHN0cmljdCc7XG5cbiFmdW5jdGlvbihleHBvcnRzLCB1bmRlZmluZWQpIHtcblxuICB2YXJcbiAgICAvLyBJZiB0aGUgdHlwZWQgYXJyYXkgaXMgdW5zcGVjaWZpZWQsIHVzZSB0aGlzLlxuICAgIERlZmF1bHRBcnJheVR5cGUgPSBGbG9hdDMyQXJyYXksXG4gICAgLy8gU2ltcGxlIG1hdGggZnVuY3Rpb25zIHdlIG5lZWQuXG4gICAgc3FydCA9IE1hdGguc3FydCxcbiAgICBzcXIgPSBmdW5jdGlvbihudW1iZXIpIHtyZXR1cm4gTWF0aC5wb3cobnVtYmVyLCAyKX0sXG4gICAgLy8gSW50ZXJuYWwgY29udmVuaWVuY2UgY29waWVzIG9mIHRoZSBleHBvcnRlZCBmdW5jdGlvbnNcbiAgICBpc0NvbXBsZXhBcnJheSxcbiAgICBDb21wbGV4QXJyYXlcblxuICBleHBvcnRzLmlzQ29tcGxleEFycmF5ID0gaXNDb21wbGV4QXJyYXkgPSBmdW5jdGlvbihvYmopIHtcbiAgICByZXR1cm4gb2JqICE9PSB1bmRlZmluZWQgJiZcbiAgICAgIG9iai5oYXNPd25Qcm9wZXJ0eSAhPT0gdW5kZWZpbmVkICYmXG4gICAgICBvYmouaGFzT3duUHJvcGVydHkoJ3JlYWwnKSAmJlxuICAgICAgb2JqLmhhc093blByb3BlcnR5KCdpbWFnJylcbiAgfVxuXG4gIGV4cG9ydHMuQ29tcGxleEFycmF5ID0gQ29tcGxleEFycmF5ID0gZnVuY3Rpb24ob3RoZXIsIG9wdF9hcnJheV90eXBlKXtcbiAgICBpZiAoaXNDb21wbGV4QXJyYXkob3RoZXIpKSB7XG4gICAgICAvLyBDb3B5IGNvbnN0dWN0b3IuXG4gICAgICB0aGlzLkFycmF5VHlwZSA9IG90aGVyLkFycmF5VHlwZVxuICAgICAgdGhpcy5yZWFsID0gbmV3IHRoaXMuQXJyYXlUeXBlKG90aGVyLnJlYWwpXG4gICAgICB0aGlzLmltYWcgPSBuZXcgdGhpcy5BcnJheVR5cGUob3RoZXIuaW1hZylcbiAgICB9IGVsc2Uge1xuICAgICAgdGhpcy5BcnJheVR5cGUgPSBvcHRfYXJyYXlfdHlwZSB8fCBEZWZhdWx0QXJyYXlUeXBlXG4gICAgICAvLyBvdGhlciBjYW4gYmUgZWl0aGVyIGFuIGFycmF5IG9yIGEgbnVtYmVyLlxuICAgICAgdGhpcy5yZWFsID0gbmV3IHRoaXMuQXJyYXlUeXBlKG90aGVyKVxuICAgICAgdGhpcy5pbWFnID0gbmV3IHRoaXMuQXJyYXlUeXBlKHRoaXMucmVhbC5sZW5ndGgpXG4gICAgfVxuXG4gICAgdGhpcy5sZW5ndGggPSB0aGlzLnJlYWwubGVuZ3RoXG4gIH1cblxuICBDb21wbGV4QXJyYXkucHJvdG90eXBlLnRvU3RyaW5nID0gZnVuY3Rpb24oKSB7XG4gICAgdmFyIGNvbXBvbmVudHMgPSBbXVxuXG4gICAgdGhpcy5mb3JFYWNoKGZ1bmN0aW9uKGNfdmFsdWUsIGkpIHtcbiAgICAgIGNvbXBvbmVudHMucHVzaChcbiAgICAgICAgJygnICtcbiAgICAgICAgY192YWx1ZS5yZWFsLnRvRml4ZWQoMikgKyAnLCcgK1xuICAgICAgICBjX3ZhbHVlLmltYWcudG9GaXhlZCgyKSArXG4gICAgICAgICcpJ1xuICAgICAgKVxuICAgIH0pXG5cbiAgICByZXR1cm4gJ1snICsgY29tcG9uZW50cy5qb2luKCcsJykgKyAnXSdcbiAgfVxuXG4gIC8vIEluLXBsYWNlIG1hcHBlci5cbiAgQ29tcGxleEFycmF5LnByb3RvdHlwZS5tYXAgPSBmdW5jdGlvbihtYXBwZXIpIHtcbiAgICB2YXJcbiAgICAgIGksXG4gICAgICBuID0gdGhpcy5sZW5ndGgsXG4gICAgICAvLyBGb3IgR0MgZWZmaWNpZW5jeSwgcGFzcyBhIHNpbmdsZSBjX3ZhbHVlIG9iamVjdCB0byB0aGUgbWFwcGVyLlxuICAgICAgY192YWx1ZSA9IHt9XG5cbiAgICBmb3IgKGkgPSAwOyBpIDwgbjsgaSsrKSB7XG4gICAgICBjX3ZhbHVlLnJlYWwgPSB0aGlzLnJlYWxbaV1cbiAgICAgIGNfdmFsdWUuaW1hZyA9IHRoaXMuaW1hZ1tpXVxuICAgICAgbWFwcGVyKGNfdmFsdWUsIGksIG4pXG4gICAgICB0aGlzLnJlYWxbaV0gPSBjX3ZhbHVlLnJlYWxcbiAgICAgIHRoaXMuaW1hZ1tpXSA9IGNfdmFsdWUuaW1hZ1xuICAgIH1cblxuICAgIHJldHVybiB0aGlzXG4gIH1cblxuICBDb21wbGV4QXJyYXkucHJvdG90eXBlLmZvckVhY2ggPSBmdW5jdGlvbihpdGVyYXRvcikge1xuICAgIHZhclxuICAgICAgaSxcbiAgICAgIG4gPSB0aGlzLmxlbmd0aCxcbiAgICAgIC8vIEZvciBjb25zaXN0ZW5jeSB3aXRoIC5tYXAuXG4gICAgICBjX3ZhbHVlID0ge31cblxuICAgIGZvciAoaSA9IDA7IGkgPCBuOyBpKyspIHtcbiAgICAgIGNfdmFsdWUucmVhbCA9IHRoaXMucmVhbFtpXVxuICAgICAgY192YWx1ZS5pbWFnID0gdGhpcy5pbWFnW2ldXG4gICAgICBpdGVyYXRvcihjX3ZhbHVlLCBpLCBuKVxuICAgIH1cbiAgfVxuXG4gIENvbXBsZXhBcnJheS5wcm90b3R5cGUuY29uanVnYXRlID0gZnVuY3Rpb24oKSB7XG4gICAgcmV0dXJuIChuZXcgQ29tcGxleEFycmF5KHRoaXMpKS5tYXAoZnVuY3Rpb24odmFsdWUpIHtcbiAgICAgIHZhbHVlLmltYWcgKj0gLTFcbiAgICB9KVxuICB9XG5cbiAgLy8gSGVscGVyIHNvIHdlIGNhbiBtYWtlIEFycmF5VHlwZSBvYmplY3RzIHJldHVybmVkIGhhdmUgc2ltaWxhciBpbnRlcmZhY2VzXG4gIC8vICAgdG8gQ29tcGxleEFycmF5cy5cbiAgZnVuY3Rpb24gaXRlcmFibGUob2JqKSB7XG4gICAgaWYgKCFvYmouZm9yRWFjaClcbiAgICAgIG9iai5mb3JFYWNoID0gZnVuY3Rpb24oaXRlcmF0b3IpIHtcbiAgICAgICAgdmFyIGksIG4gPSB0aGlzLmxlbmd0aFxuXG4gICAgICAgIGZvciAoaSA9IDA7IGkgPCBuOyBpKyspXG4gICAgICAgICAgaXRlcmF0b3IodGhpc1tpXSwgaSwgbilcbiAgICAgIH1cblxuICAgIHJldHVybiBvYmpcbiAgfVxuXG4gIENvbXBsZXhBcnJheS5wcm90b3R5cGUubWFnbml0dWRlID0gZnVuY3Rpb24oKSB7XG4gICAgdmFyIG1hZ3MgPSBuZXcgdGhpcy5BcnJheVR5cGUodGhpcy5sZW5ndGgpXG5cbiAgICB0aGlzLmZvckVhY2goZnVuY3Rpb24odmFsdWUsIGkpIHtcbiAgICAgIG1hZ3NbaV0gPSBzcXJ0KHNxcih2YWx1ZS5yZWFsKSArIHNxcih2YWx1ZS5pbWFnKSlcbiAgICB9KVxuXG4gICAgLy8gQXJyYXlUeXBlIHdpbGwgbm90IG5lY2Vzc2FyaWx5IGJlIGl0ZXJhYmxlOiBtYWtlIGl0IHNvLlxuICAgIHJldHVybiBpdGVyYWJsZShtYWdzKVxuICB9XG59KHR5cGVvZiBleHBvcnRzID09PSAndW5kZWZpbmVkJyAmJiAodGhpcy5jb21wbGV4X2FycmF5ID0ge30pIHx8IGV4cG9ydHMpXG4iLCIndXNlIHN0cmljdCc7XG5cbiFmdW5jdGlvbihleHBvcnRzLCBjb21wbGV4X2FycmF5KSB7XG5cbiAgdmFyXG4gICAgQ29tcGxleEFycmF5ID0gY29tcGxleF9hcnJheS5Db21wbGV4QXJyYXksXG4gICAgLy8gTWF0aCBjb25zdGFudHMgYW5kIGZ1bmN0aW9ucyB3ZSBuZWVkLlxuICAgIFBJID0gTWF0aC5QSSxcbiAgICBTUVJUMV8yID0gTWF0aC5TUVJUMV8yLFxuICAgIHNxcnQgPSBNYXRoLnNxcnQsXG4gICAgY29zID0gTWF0aC5jb3MsXG4gICAgc2luID0gTWF0aC5zaW5cblxuICBDb21wbGV4QXJyYXkucHJvdG90eXBlLkZGVCA9IGZ1bmN0aW9uKCkge1xuICAgIHJldHVybiBGRlQodGhpcywgZmFsc2UpO1xuICB9XG5cbiAgZXhwb3J0cy5GRlQgPSBmdW5jdGlvbihpbnB1dCkge1xuICAgIHJldHVybiBlbnN1cmVDb21wbGV4QXJyYXkoaW5wdXQpLkZGVCgpXG4gIH1cblxuICBDb21wbGV4QXJyYXkucHJvdG90eXBlLkludkZGVCA9IGZ1bmN0aW9uKCkge1xuICAgIHJldHVybiBGRlQodGhpcywgdHJ1ZSlcbiAgfVxuXG4gIGV4cG9ydHMuSW52RkZUID0gZnVuY3Rpb24oaW5wdXQpIHtcbiAgICByZXR1cm4gZW5zdXJlQ29tcGxleEFycmF5KGlucHV0KS5JbnZGRlQoKVxuICB9XG5cbiAgLy8gQXBwbGllcyBhIGZyZXF1ZW5jeS1zcGFjZSBmaWx0ZXIgdG8gaW5wdXQsIGFuZCByZXR1cm5zIHRoZSByZWFsLXNwYWNlXG4gIC8vIGZpbHRlcmVkIGlucHV0LlxuICAvLyBmaWx0ZXJlciBhY2NlcHRzIGZyZXEsIGksIG4gYW5kIG1vZGlmaWVzIGZyZXEucmVhbCBhbmQgZnJlcS5pbWFnLlxuICBDb21wbGV4QXJyYXkucHJvdG90eXBlLmZyZXF1ZW5jeU1hcCA9IGZ1bmN0aW9uKGZpbHRlcmVyKSB7XG4gICAgcmV0dXJuIHRoaXMuRkZUKCkubWFwKGZpbHRlcmVyKS5JbnZGRlQoKVxuICB9XG5cbiAgZXhwb3J0cy5mcmVxdWVuY3lNYXAgPSBmdW5jdGlvbihpbnB1dCwgZmlsdGVyZXIpIHtcbiAgICByZXR1cm4gZW5zdXJlQ29tcGxleEFycmF5KGlucHV0KS5mcmVxdWVuY3lNYXAoZmlsdGVyZXIpXG4gIH1cblxuICBmdW5jdGlvbiBlbnN1cmVDb21wbGV4QXJyYXkoaW5wdXQpIHtcbiAgICByZXR1cm4gY29tcGxleF9hcnJheS5pc0NvbXBsZXhBcnJheShpbnB1dCkgJiYgaW5wdXQgfHxcbiAgICAgICAgbmV3IENvbXBsZXhBcnJheShpbnB1dClcbiAgfVxuXG4gIGZ1bmN0aW9uIEZGVChpbnB1dCwgaW52ZXJzZSkge1xuICAgIHZhciBuID0gaW5wdXQubGVuZ3RoXG5cbiAgICBpZiAobiAmIChuIC0gMSkpIHtcbiAgICAgIHJldHVybiBGRlRfUmVjdXJzaXZlKGlucHV0LCBpbnZlcnNlKVxuICAgIH0gZWxzZSB7XG4gICAgICByZXR1cm4gRkZUXzJfSXRlcmF0aXZlKGlucHV0LCBpbnZlcnNlKVxuICAgIH1cbiAgfVxuXG4gIGZ1bmN0aW9uIEZGVF9SZWN1cnNpdmUoaW5wdXQsIGludmVyc2UpIHtcbiAgICB2YXJcbiAgICAgIG4gPSBpbnB1dC5sZW5ndGgsXG4gICAgICAvLyBDb3VudGVycy5cbiAgICAgIGksIGosXG4gICAgICBvdXRwdXQsXG4gICAgICAvLyBDb21wbGV4IG11bHRpcGxpZXIgYW5kIGl0cyBkZWx0YS5cbiAgICAgIGZfciwgZl9pLCBkZWxfZl9yLCBkZWxfZl9pLFxuICAgICAgLy8gTG93ZXN0IGRpdmlzb3IgYW5kIHJlbWFpbmRlci5cbiAgICAgIHAsIG0sXG4gICAgICBub3JtYWxpc2F0aW9uLFxuICAgICAgcmVjdXJzaXZlX3Jlc3VsdCxcbiAgICAgIF9zd2FwLCBfcmVhbCwgX2ltYWdcblxuICAgIGlmIChuID09PSAxKSB7XG4gICAgICByZXR1cm4gaW5wdXRcbiAgICB9XG5cbiAgICBvdXRwdXQgPSBuZXcgQ29tcGxleEFycmF5KG4sIGlucHV0LkFycmF5VHlwZSlcblxuICAgIC8vIFVzZSB0aGUgbG93ZXN0IG9kZCBmYWN0b3IsIHNvIHdlIGFyZSBhYmxlIHRvIHVzZSBGRlRfMl9JdGVyYXRpdmUgaW4gdGhlXG4gICAgLy8gcmVjdXJzaXZlIHRyYW5zZm9ybXMgb3B0aW1hbGx5LlxuICAgIHAgPSBMb3dlc3RPZGRGYWN0b3IobilcbiAgICBtID0gbiAvIHBcbiAgICBub3JtYWxpc2F0aW9uID0gMSAvIHNxcnQocClcbiAgICByZWN1cnNpdmVfcmVzdWx0ID0gbmV3IENvbXBsZXhBcnJheShtLCBpbnB1dC5BcnJheVR5cGUpXG5cbiAgICAvLyBMb29wcyBnbyBsaWtlIE8obiDOoyBwX2kpLCB3aGVyZSBwX2kgYXJlIHRoZSBwcmltZSBmYWN0b3JzIG9mIG4uXG4gICAgLy8gZm9yIGEgcG93ZXIgb2YgYSBwcmltZSwgcCwgdGhpcyByZWR1Y2VzIHRvIE8obiBwIGxvZ19wIG4pXG4gICAgZm9yKGogPSAwOyBqIDwgcDsgaisrKSB7XG4gICAgICBmb3IoaSA9IDA7IGkgPCBtOyBpKyspIHtcbiAgICAgICAgcmVjdXJzaXZlX3Jlc3VsdC5yZWFsW2ldID0gaW5wdXQucmVhbFtpICogcCArIGpdXG4gICAgICAgIHJlY3Vyc2l2ZV9yZXN1bHQuaW1hZ1tpXSA9IGlucHV0LmltYWdbaSAqIHAgKyBqXVxuICAgICAgfVxuICAgICAgLy8gRG9uJ3QgZ28gZGVlcGVyIHVubGVzcyBuZWNlc3NhcnkgdG8gc2F2ZSBhbGxvY3MuXG4gICAgICBpZiAobSA+IDEpIHtcbiAgICAgICAgcmVjdXJzaXZlX3Jlc3VsdCA9IEZGVChyZWN1cnNpdmVfcmVzdWx0LCBpbnZlcnNlKVxuICAgICAgfVxuXG4gICAgICBkZWxfZl9yID0gY29zKDIqUEkqai9uKVxuICAgICAgZGVsX2ZfaSA9IChpbnZlcnNlID8gLTEgOiAxKSAqIHNpbigyKlBJKmovbilcbiAgICAgIGZfciA9IDFcbiAgICAgIGZfaSA9IDBcblxuICAgICAgZm9yKGkgPSAwOyBpIDwgbjsgaSsrKSB7XG4gICAgICAgIF9yZWFsID0gcmVjdXJzaXZlX3Jlc3VsdC5yZWFsW2kgJSBtXVxuICAgICAgICBfaW1hZyA9IHJlY3Vyc2l2ZV9yZXN1bHQuaW1hZ1tpICUgbV1cblxuICAgICAgICBvdXRwdXQucmVhbFtpXSArPSBmX3IgKiBfcmVhbCAtIGZfaSAqIF9pbWFnXG4gICAgICAgIG91dHB1dC5pbWFnW2ldICs9IGZfciAqIF9pbWFnICsgZl9pICogX3JlYWxcblxuICAgICAgICBfc3dhcCA9IGZfciAqIGRlbF9mX3IgLSBmX2kgKiBkZWxfZl9pXG4gICAgICAgIGZfaSA9IGZfciAqIGRlbF9mX2kgKyBmX2kgKiBkZWxfZl9yXG4gICAgICAgIGZfciA9IF9zd2FwXG4gICAgICB9XG4gICAgfVxuXG4gICAgLy8gQ29weSBiYWNrIHRvIGlucHV0IHRvIG1hdGNoIEZGVF8yX0l0ZXJhdGl2ZSBpbi1wbGFjZW5lc3NcbiAgICAvLyBUT0RPOiBmYXN0ZXIgd2F5IG9mIG1ha2luZyB0aGlzIGluLXBsYWNlP1xuICAgIGZvcihpID0gMDsgaSA8IG47IGkrKykge1xuICAgICAgaW5wdXQucmVhbFtpXSA9IG5vcm1hbGlzYXRpb24gKiBvdXRwdXQucmVhbFtpXVxuICAgICAgaW5wdXQuaW1hZ1tpXSA9IG5vcm1hbGlzYXRpb24gKiBvdXRwdXQuaW1hZ1tpXVxuICAgIH1cblxuICAgIHJldHVybiBpbnB1dFxuICB9XG5cbiAgZnVuY3Rpb24gRkZUXzJfSXRlcmF0aXZlKGlucHV0LCBpbnZlcnNlKSB7XG4gICAgdmFyXG4gICAgICBuID0gaW5wdXQubGVuZ3RoLFxuICAgICAgLy8gQ291bnRlcnMuXG4gICAgICBpLCBqLFxuICAgICAgb3V0cHV0LCBvdXRwdXRfciwgb3V0cHV0X2ksXG4gICAgICAvLyBDb21wbGV4IG11bHRpcGxpZXIgYW5kIGl0cyBkZWx0YS5cbiAgICAgIGZfciwgZl9pLCBkZWxfZl9yLCBkZWxfZl9pLCB0ZW1wLFxuICAgICAgLy8gVGVtcG9yYXJ5IGxvb3AgdmFyaWFibGVzLlxuICAgICAgbF9pbmRleCwgcl9pbmRleCxcbiAgICAgIGxlZnRfciwgbGVmdF9pLCByaWdodF9yLCByaWdodF9pLFxuICAgICAgLy8gd2lkdGggb2YgZWFjaCBzdWItYXJyYXkgZm9yIHdoaWNoIHdlJ3JlIGl0ZXJhdGl2ZWx5IGNhbGN1bGF0aW5nIEZGVC5cbiAgICAgIHdpZHRoXG5cbiAgICBvdXRwdXQgPSBCaXRSZXZlcnNlQ29tcGxleEFycmF5KGlucHV0KVxuICAgIG91dHB1dF9yID0gb3V0cHV0LnJlYWxcbiAgICBvdXRwdXRfaSA9IG91dHB1dC5pbWFnXG4gICAgLy8gTG9vcHMgZ28gbGlrZSBPKG4gbG9nIG4pOlxuICAgIC8vICAgd2lkdGggfiBsb2cgbjsgaSxqIH4gblxuICAgIHdpZHRoID0gMVxuICAgIHdoaWxlICh3aWR0aCA8IG4pIHtcbiAgICAgIGRlbF9mX3IgPSBjb3MoUEkvd2lkdGgpXG4gICAgICBkZWxfZl9pID0gKGludmVyc2UgPyAtMSA6IDEpICogc2luKFBJL3dpZHRoKVxuICAgICAgZm9yIChpID0gMDsgaSA8IG4vKDIqd2lkdGgpOyBpKyspIHtcbiAgICAgICAgZl9yID0gMVxuICAgICAgICBmX2kgPSAwXG4gICAgICAgIGZvciAoaiA9IDA7IGogPCB3aWR0aDsgaisrKSB7XG4gICAgICAgICAgbF9pbmRleCA9IDIqaSp3aWR0aCArIGpcbiAgICAgICAgICByX2luZGV4ID0gbF9pbmRleCArIHdpZHRoXG5cbiAgICAgICAgICBsZWZ0X3IgPSBvdXRwdXRfcltsX2luZGV4XVxuICAgICAgICAgIGxlZnRfaSA9IG91dHB1dF9pW2xfaW5kZXhdXG4gICAgICAgICAgcmlnaHRfciA9IGZfciAqIG91dHB1dF9yW3JfaW5kZXhdIC0gZl9pICogb3V0cHV0X2lbcl9pbmRleF1cbiAgICAgICAgICByaWdodF9pID0gZl9pICogb3V0cHV0X3Jbcl9pbmRleF0gKyBmX3IgKiBvdXRwdXRfaVtyX2luZGV4XVxuXG4gICAgICAgICAgb3V0cHV0X3JbbF9pbmRleF0gPSBTUVJUMV8yICogKGxlZnRfciArIHJpZ2h0X3IpXG4gICAgICAgICAgb3V0cHV0X2lbbF9pbmRleF0gPSBTUVJUMV8yICogKGxlZnRfaSArIHJpZ2h0X2kpXG4gICAgICAgICAgb3V0cHV0X3Jbcl9pbmRleF0gPSBTUVJUMV8yICogKGxlZnRfciAtIHJpZ2h0X3IpXG4gICAgICAgICAgb3V0cHV0X2lbcl9pbmRleF0gPSBTUVJUMV8yICogKGxlZnRfaSAtIHJpZ2h0X2kpXG4gICAgICAgICAgdGVtcCA9IGZfciAqIGRlbF9mX3IgLSBmX2kgKiBkZWxfZl9pXG4gICAgICAgICAgZl9pID0gZl9yICogZGVsX2ZfaSArIGZfaSAqIGRlbF9mX3JcbiAgICAgICAgICBmX3IgPSB0ZW1wXG4gICAgICAgIH1cbiAgICAgIH1cbiAgICAgIHdpZHRoIDw8PSAxXG4gICAgfVxuXG4gICAgcmV0dXJuIG91dHB1dFxuICB9XG5cbiAgZnVuY3Rpb24gQml0UmV2ZXJzZUluZGV4KGluZGV4LCBuKSB7XG4gICAgdmFyIGJpdHJldmVyc2VkX2luZGV4ID0gMFxuXG4gICAgd2hpbGUgKG4gPiAxKSB7XG4gICAgICBiaXRyZXZlcnNlZF9pbmRleCA8PD0gMVxuICAgICAgYml0cmV2ZXJzZWRfaW5kZXggKz0gaW5kZXggJiAxXG4gICAgICBpbmRleCA+Pj0gMVxuICAgICAgbiA+Pj0gMVxuICAgIH1cbiAgICByZXR1cm4gYml0cmV2ZXJzZWRfaW5kZXhcbiAgfVxuXG4gIGZ1bmN0aW9uIEJpdFJldmVyc2VDb21wbGV4QXJyYXkoYXJyYXkpIHtcbiAgICB2YXIgbiA9IGFycmF5Lmxlbmd0aCxcbiAgICAgICAgZmxpcHMgPSB7fSxcbiAgICAgICAgc3dhcCxcbiAgICAgICAgaVxuXG4gICAgZm9yKGkgPSAwOyBpIDwgbjsgaSsrKSB7XG4gICAgICB2YXIgcl9pID0gQml0UmV2ZXJzZUluZGV4KGksIG4pXG5cbiAgICAgIGlmIChmbGlwcy5oYXNPd25Qcm9wZXJ0eShpKSB8fCBmbGlwcy5oYXNPd25Qcm9wZXJ0eShyX2kpKSBjb250aW51ZVxuXG4gICAgICBzd2FwID0gYXJyYXkucmVhbFtyX2ldXG4gICAgICBhcnJheS5yZWFsW3JfaV0gPSBhcnJheS5yZWFsW2ldXG4gICAgICBhcnJheS5yZWFsW2ldID0gc3dhcFxuXG4gICAgICBzd2FwID0gYXJyYXkuaW1hZ1tyX2ldXG4gICAgICBhcnJheS5pbWFnW3JfaV0gPSBhcnJheS5pbWFnW2ldXG4gICAgICBhcnJheS5pbWFnW2ldID0gc3dhcFxuXG4gICAgICBmbGlwc1tpXSA9IGZsaXBzW3JfaV0gPSB0cnVlXG4gICAgfVxuXG4gICAgcmV0dXJuIGFycmF5XG4gIH1cblxuICBmdW5jdGlvbiBMb3dlc3RPZGRGYWN0b3Iobikge1xuICAgIHZhciBmYWN0b3IgPSAzLFxuICAgICAgICBzcXJ0X24gPSBzcXJ0KG4pXG5cbiAgICB3aGlsZShmYWN0b3IgPD0gc3FydF9uKSB7XG4gICAgICBpZiAobiAlIGZhY3RvciA9PT0gMCkgcmV0dXJuIGZhY3RvclxuICAgICAgZmFjdG9yID0gZmFjdG9yICsgMlxuICAgIH1cbiAgICByZXR1cm4gblxuICB9XG5cbn0oXG4gIHR5cGVvZiBleHBvcnRzID09PSAndW5kZWZpbmVkJyAmJiAodGhpcy5mZnQgPSB7fSkgfHwgZXhwb3J0cyxcbiAgdHlwZW9mIHJlcXVpcmUgPT09ICd1bmRlZmluZWQnICYmICh0aGlzLmNvbXBsZXhfYXJyYXkpIHx8XG4gICAgcmVxdWlyZSgnLi9jb21wbGV4X2FycmF5JylcbilcbiIsIm1vZHVsZS5leHBvcnRzID0gZnVuY3Rpb24oYnVmZmVyU2l6ZSwgbSl7XG4gIHJldHVybiBtLmFtcFNwZWN0cnVtO1xufTsiLCJtb2R1bGUuZXhwb3J0cyA9IGZ1bmN0aW9uKGJ1ZmZlclNpemUsIG0pIHtcbiAgcmV0dXJuIG0uc2lnbmFsO1xufTsiLCJtb2R1bGUuZXhwb3J0cyA9IGZ1bmN0aW9uKGJ1ZmZlclNpemUsIG0pIHtcbiAgcmV0dXJuIG0uY29tcGxleFNwZWN0cnVtO1xufTsiLCJtb2R1bGUuZXhwb3J0cyA9IGZ1bmN0aW9uKGJ1ZmZlclNpemUsIG0pIHtcbiAgdmFyIGVuZXJneSA9IDA7XG4gIGZvciAodmFyIGkgPSAwOyBpIDwgbS5zaWduYWwubGVuZ3RoOyBpKyspIHtcbiAgICBlbmVyZ3kgKz0gTWF0aC5wb3coTWF0aC5hYnMobS5zaWduYWxbaV0pLCAyKTtcbiAgfVxuICByZXR1cm4gZW5lcmd5O1xufTsiLCJtb2R1bGUuZXhwb3J0cyA9IHtcbiAgYnVmZmVyOiByZXF1aXJlKCcuL2J1ZmZlcicpLFxuICBybXM6IHJlcXVpcmUoJy4vcm1zJyksXG4gIGVuZXJneTogcmVxdWlyZSgnLi9lbmVyZ3knKSxcbiAgY29tcGxleFNwZWN0cnVtOiByZXF1aXJlKCcuL2NvbXBsZXhTcGVjdHJ1bScpLFxuICBzcGVjdHJhbFNsb3BlOiByZXF1aXJlKCcuL3NwZWN0cmFsU2xvcGUnKSxcbiAgc3BlY3RyYWxDZW50cm9pZDogcmVxdWlyZSgnLi9zcGVjdHJhbENlbnRyb2lkJyksXG4gIHNwZWN0cmFsUm9sbG9mZjogcmVxdWlyZSgnLi9zcGVjdHJhbFJvbGxvZmYnKSxcbiAgc3BlY3RyYWxGbGF0bmVzczogcmVxdWlyZSgnLi9zcGVjdHJhbEZsYXRuZXNzJyksXG4gIHNwZWN0cmFsU3ByZWFkOiByZXF1aXJlKCcuL3NwZWN0cmFsU3ByZWFkJyksXG4gIHNwZWN0cmFsU2tld25lc3M6IHJlcXVpcmUoJy4vc3BlY3RyYWxTa2V3bmVzcycpLFxuICBzcGVjdHJhbEt1cnRvc2lzOiByZXF1aXJlKCcuL3NwZWN0cmFsS3VydG9zaXMnKSxcbiAgYW1wbGl0dWRlU3BlY3RydW06IHJlcXVpcmUoJy4vYW1wbGl0dWRlU3BlY3RydW0nKSxcbiAgemNyOiByZXF1aXJlKCcuL3pjcicpLFxuICBwb3dlclNwZWN0cnVtOiByZXF1aXJlKCcuL3Bvd2VyU3BlY3RydW0nKSxcbiAgbG91ZG5lc3M6IHJlcXVpcmUoJy4vbG91ZG5lc3MnKSxcbiAgcGVyY2VwdHVhbFNwcmVhZDogcmVxdWlyZSgnLi9wZXJjZXB0dWFsU3ByZWFkJyksXG4gIHBlcmNlcHR1YWxTaGFycG5lc3M6IHJlcXVpcmUoJy4vcGVyY2VwdHVhbFNoYXJwbmVzcycpLFxuICBtZmNjOiByZXF1aXJlKCcuL21mY2MnKVxufTsiLCJtb2R1bGUuZXhwb3J0cyA9IGZ1bmN0aW9uKGJ1ZmZlclNpemUsIG0pIHtcbiAgdmFyIGJhcmtTY2FsZSA9IG5ldyBGbG9hdDMyQXJyYXkobS5hbXBTcGVjdHJ1bS5sZW5ndGgpO1xuICB2YXIgTlVNX0JBUktfQkFORFMgPSAyNDtcbiAgdmFyIHNwZWNpZmljID0gbmV3IEZsb2F0MzJBcnJheShOVU1fQkFSS19CQU5EUyk7XG4gIHZhciB0b3QgPSAwO1xuICB2YXIgbm9ybWFsaXNlZFNwZWN0cnVtID0gbS5hbXBTcGVjdHJ1bTtcbiAgdmFyIGJiTGltaXRzID0gbmV3IEludDMyQXJyYXkoTlVNX0JBUktfQkFORFMgKyAxKTtcblxuICBmb3IgKHZhciBpID0gMDsgaSA8IGJhcmtTY2FsZS5sZW5ndGg7IGkrKykge1xuICAgIGJhcmtTY2FsZVtpXSA9IGkgKiBtLmF1ZGlvQ29udGV4dC5zYW1wbGVSYXRlIC8gKGJ1ZmZlclNpemUpO1xuICAgIGJhcmtTY2FsZVtpXSA9IDEzICogTWF0aC5hdGFuKGJhcmtTY2FsZVtpXSAvIDEzMTUuOCkgKyAzLjUgKiBNYXRoLmF0YW4oTWF0aC5wb3coKGJhcmtTY2FsZVtpXSAvIDc1MTgpLCAyKSk7XG4gIH1cblxuXG4gIGJiTGltaXRzWzBdID0gMDtcbiAgdmFyIGN1cnJlbnRCYW5kRW5kID0gYmFya1NjYWxlW20uYW1wU3BlY3RydW0ubGVuZ3RoIC0gMV0gLyBOVU1fQkFSS19CQU5EUztcbiAgdmFyIGN1cnJlbnRCYW5kID0gMTtcbiAgZm9yICh2YXIgaSA9IDA7IGkgPCBtLmFtcFNwZWN0cnVtLmxlbmd0aDsgaSsrKSB7XG4gICAgd2hpbGUgKGJhcmtTY2FsZVtpXSA+IGN1cnJlbnRCYW5kRW5kKSB7XG4gICAgICBiYkxpbWl0c1tjdXJyZW50QmFuZCsrXSA9IGk7XG4gICAgICBjdXJyZW50QmFuZEVuZCA9IGN1cnJlbnRCYW5kICogYmFya1NjYWxlW20uYW1wU3BlY3RydW0ubGVuZ3RoIC0gMV0gLyBOVU1fQkFSS19CQU5EUztcbiAgICB9XG4gIH1cblxuICBiYkxpbWl0c1tOVU1fQkFSS19CQU5EU10gPSBtLmFtcFNwZWN0cnVtLmxlbmd0aCAtIDE7XG5cbiAgLy9wcm9jZXNzXG5cbiAgZm9yICh2YXIgaSA9IDA7IGkgPCBOVU1fQkFSS19CQU5EUzsgaSsrKSB7XG4gICAgdmFyIHN1bSA9IDA7XG4gICAgZm9yICh2YXIgaiA9IGJiTGltaXRzW2ldOyBqIDwgYmJMaW1pdHNbaSArIDFdOyBqKyspIHtcblxuICAgICAgc3VtICs9IG5vcm1hbGlzZWRTcGVjdHJ1bVtqXTtcbiAgICB9XG4gICAgc3BlY2lmaWNbaV0gPSBNYXRoLnBvdyhzdW0sIDAuMjMpO1xuICB9XG5cbiAgLy9nZXQgdG90YWwgbG91ZG5lc3NcbiAgZm9yICh2YXIgaSA9IDA7IGkgPCBzcGVjaWZpYy5sZW5ndGg7IGkrKykge1xuICAgIHRvdCArPSBzcGVjaWZpY1tpXTtcbiAgfVxuICByZXR1cm4ge1xuICAgIFwic3BlY2lmaWNcIjogc3BlY2lmaWMsXG4gICAgXCJ0b3RhbFwiOiB0b3RcbiAgfTtcbn07IiwiLy91c2VkIHR1dG9yaWFsIGZyb20gaHR0cDovL3ByYWN0aWNhbGNyeXB0b2dyYXBoeS5jb20vbWlzY2VsbGFuZW91cy9tYWNoaW5lLWxlYXJuaW5nL2d1aWRlLW1lbC1mcmVxdWVuY3ktY2Vwc3RyYWwtY29lZmZpY2llbnRzLW1mY2NzL1xuXG52YXIgcG93ZXJTcGVjdHJ1bSA9IHJlcXVpcmUoJy4vcG93ZXJTcGVjdHJ1bScpO1xuXG5tb2R1bGUuZXhwb3J0cyA9IGZ1bmN0aW9uKGJ1ZmZlclNpemUsIG0pIHtcbiAgdmFyIHBvd1NwZWMgPSBwb3dlclNwZWN0cnVtKGJ1ZmZlclNpemUsIG0pO1xuICB2YXIgZnJlcVRvTWVsID0gZnVuY3Rpb24oZnJlcVZhbHVlKSB7XG4gICAgdmFyIG1lbFZhbHVlID0gMTEyNSAqIE1hdGgubG9nKDEgKyAoZnJlcVZhbHVlIC8gNzAwKSk7XG4gICAgcmV0dXJuIG1lbFZhbHVlO1xuICB9O1xuICB2YXIgbWVsVG9GcmVxID0gZnVuY3Rpb24obWVsVmFsdWUpIHtcbiAgICB2YXIgZnJlcVZhbHVlID0gNzAwICogKE1hdGguZXhwKG1lbFZhbHVlIC8gMTEyNSkgLSAxKTtcbiAgICByZXR1cm4gZnJlcVZhbHVlO1xuICB9O1xuICB2YXIgbnVtRmlsdGVycyA9IDI2OyAvLzI2IGZpbHRlcnMgaXMgc3RhbmRhcmRcbiAgdmFyIG1lbFZhbHVlcyA9IG5ldyBGbG9hdDMyQXJyYXkobnVtRmlsdGVycyArIDIpOyAvL3RoZSArMiBpcyB0aGUgdXBwZXIgYW5kIGxvd2VyIGxpbWl0c1xuICB2YXIgbWVsVmFsdWVzSW5GcmVxID0gbmV3IEZsb2F0MzJBcnJheShudW1GaWx0ZXJzICsgMik7XG4gIC8vR2VuZXJhdGUgbGltaXRzIGluIEh6IC0gZnJvbSAwIHRvIHRoZSBueXF1aXN0LlxuICB2YXIgbG93ZXJMaW1pdEZyZXEgPSAwO1xuICB2YXIgdXBwZXJMaW1pdEZyZXEgPSBhdWRpb0NvbnRleHQuc2FtcGxlUmF0ZSAvIDI7XG4gIC8vQ29udmVydCB0aGUgbGltaXRzIHRvIE1lbFxuICB2YXIgbG93ZXJMaW1pdE1lbCA9IGZyZXFUb01lbChsb3dlckxpbWl0RnJlcSk7XG4gIHZhciB1cHBlckxpbWl0TWVsID0gZnJlcVRvTWVsKHVwcGVyTGltaXRGcmVxKTtcbiAgLy9GaW5kIHRoZSByYW5nZVxuICB2YXIgcmFuZ2UgPSB1cHBlckxpbWl0TWVsIC0gbG93ZXJMaW1pdE1lbDtcbiAgLy9GaW5kIHRoZSByYW5nZSBhcyBwYXJ0IG9mIHRoZSBsaW5lYXIgaW50ZXJwb2xhdGlvblxuICB2YXIgdmFsdWVUb0FkZCA9IHJhbmdlIC8gKG51bUZpbHRlcnMgKyAxKTtcblxuICB2YXIgZmZ0Qmluc09mRnJlcSA9IEFycmF5KG51bUZpbHRlcnMgKyAyKTtcblxuICBmb3IgKHZhciBpID0gMDsgaSA8IG1lbFZhbHVlcy5sZW5ndGg7IGkrKykge1xuICAgIC8vSW5pdGlhbGlzaW5nIHRoZSBtZWwgZnJlcXVlbmNpZXMgLSB0aGV5IGFyZSBqdXN0IGEgbGluZWFyIGludGVycG9sYXRpb24gYmV0d2VlbiB0aGUgbG93ZXIgYW5kIHVwcGVyIGxpbWl0cy5cbiAgICBtZWxWYWx1ZXNbaV0gPSBpICogdmFsdWVUb0FkZDtcbiAgICAvL0NvbnZlcnQgYmFjayB0byBIelxuICAgIG1lbFZhbHVlc0luRnJlcVtpXSA9IG1lbFRvRnJlcShtZWxWYWx1ZXNbaV0pO1xuICAgIC8vRmluZCB0aGUgY29ycmVzcG9uZGluZyBiaW5zXG4gICAgZmZ0Qmluc09mRnJlcVtpXSA9IE1hdGguZmxvb3IoKGJ1ZmZlclNpemUgKyAxKSAqIG1lbFZhbHVlc0luRnJlcVtpXSAvIGF1ZGlvQ29udGV4dC5zYW1wbGVSYXRlKTtcbiAgfVxuXG4gIHZhciBmaWx0ZXJCYW5rID0gQXJyYXkobnVtRmlsdGVycyk7XG4gIGZvciAodmFyIGogPSAwOyBqIDwgZmlsdGVyQmFuay5sZW5ndGg7IGorKykge1xuICAgIC8vY3JlYXRpbmcgYSB0d28gZGltZW5zaW9uYWwgYXJyYXkgb2Ygc2l6ZSBudW1GaWx0ZXMgKiAoYnVmZmVyc2l6ZS8yKSsxIGFuZCBwcmUtcG9wdWxhdGluZyB0aGUgYXJyYXlzIHdpdGggMHMuXG4gICAgZmlsdGVyQmFua1tqXSA9IEFycmF5LmFwcGx5KG51bGwsIG5ldyBBcnJheSgoYnVmZmVyU2l6ZSAvIDIpICsgMSkpLm1hcChOdW1iZXIucHJvdG90eXBlLnZhbHVlT2YsIDApO1xuICAgIC8vY3JlYXRpbmcgdGhlIGxvd2VyIGFuZCB1cHBlciBzbG9wZXMgZm9yIGVhY2ggYmluXG4gICAgZm9yICh2YXIgaSA9IGZmdEJpbnNPZkZyZXFbal07IGkgPCBmZnRCaW5zT2ZGcmVxW2ogKyAxXTsgaSsrKSB7XG4gICAgICBmaWx0ZXJCYW5rW2pdW2ldID0gKGkgLSBmZnRCaW5zT2ZGcmVxW2pdKSAvIChmZnRCaW5zT2ZGcmVxW2ogKyAxXSAtIGZmdEJpbnNPZkZyZXFbal0pO1xuICAgIH1cbiAgICBmb3IgKHZhciBpID0gZmZ0Qmluc09mRnJlcVtqICsgMV07IGkgPCBmZnRCaW5zT2ZGcmVxW2ogKyAyXTsgaSsrKSB7XG4gICAgICBmaWx0ZXJCYW5rW2pdW2ldID0gKGZmdEJpbnNPZkZyZXFbaiArIDJdIC0gaSkgLyAoZmZ0Qmluc09mRnJlcVtqICsgMl0gLSBmZnRCaW5zT2ZGcmVxW2ogKyAxXSk7XG4gICAgfVxuICB9XG5cbiAgdmFyIGxvZ2dlZE1lbEJhbmRzID0gbmV3IEZsb2F0MzJBcnJheShudW1GaWx0ZXJzKTtcbiAgZm9yICh2YXIgaSA9IDA7IGkgPCBsb2dnZWRNZWxCYW5kcy5sZW5ndGg7IGkrKykge1xuICAgIGxvZ2dlZE1lbEJhbmRzW2ldID0gMDtcbiAgICBmb3IgKHZhciBqID0gMDsgaiA8IChidWZmZXJTaXplIC8gMik7IGorKykge1xuICAgICAgLy9wb2ludCBtdWx0aXBsaWNhdGlvbiBiZXR3ZWVuIHBvd2VyIHNwZWN0cnVtIGFuZCBmaWx0ZXJiYW5rcy5cbiAgICAgIGZpbHRlckJhbmtbaV1bal0gPSBmaWx0ZXJCYW5rW2ldW2pdICogcG93U3BlY1tqXTtcblxuICAgICAgLy9zdW1taW5nIHVwIGFsbCBvZiB0aGUgY29lZmZpY2llbnRzIGludG8gb25lIGFycmF5XG4gICAgICBsb2dnZWRNZWxCYW5kc1tpXSArPSBmaWx0ZXJCYW5rW2ldW2pdO1xuICAgIH1cbiAgICAvL2xvZyBlYWNoIGNvZWZmaWNpZW50XG4gICAgbG9nZ2VkTWVsQmFuZHNbaV0gPSBNYXRoLmxvZyhsb2dnZWRNZWxCYW5kc1tpXSk7XG4gIH1cblxuICAvL2RjdFxuICB2YXIgayA9IE1hdGguUEkgLyBudW1GaWx0ZXJzO1xuICB2YXIgdzEgPSAxLjAgLyBNYXRoLnNxcnQobnVtRmlsdGVycyk7XG4gIHZhciB3MiA9IE1hdGguc3FydCgyLjAgLyBudW1GaWx0ZXJzKTtcbiAgdmFyIG51bUNvZWZmcyA9IDEzO1xuICB2YXIgZGN0TWF0cml4ID0gbmV3IEZsb2F0MzJBcnJheShudW1Db2VmZnMgKiBudW1GaWx0ZXJzKTtcblxuICBmb3IgKHZhciBpID0gMDsgaSA8IG51bUNvZWZmczsgaSsrKSB7XG4gICAgZm9yICh2YXIgaiA9IDA7IGogPCBudW1GaWx0ZXJzOyBqKyspIHtcbiAgICAgIHZhciBpZHggPSBpICsgKGogKiBudW1Db2VmZnMpO1xuICAgICAgaWYgKGkgPT0gMCkge1xuICAgICAgICBkY3RNYXRyaXhbaWR4XSA9IHcxICogTWF0aC5jb3MoayAqIChpICsgMSkgKiAoaiArIDAuNSkpO1xuICAgICAgfSBlbHNlIHtcbiAgICAgICAgZGN0TWF0cml4W2lkeF0gPSB3MiAqIE1hdGguY29zKGsgKiAoaSArIDEpICogKGogKyAwLjUpKTtcbiAgICAgIH1cbiAgICB9XG4gIH1cblxuICB2YXIgbWZjY3MgPSBuZXcgRmxvYXQzMkFycmF5KG51bUNvZWZmcyk7XG4gIGZvciAodmFyIGsgPSAwOyBrIDwgbnVtQ29lZmZzOyBrKyspIHtcbiAgICB2YXIgdiA9IDA7XG4gICAgZm9yICh2YXIgbiA9IDA7IG4gPCBudW1GaWx0ZXJzOyBuKyspIHtcbiAgICAgIHZhciBpZHggPSBrICsgKG4gKiBudW1Db2VmZnMpO1xuICAgICAgdiArPSAoZGN0TWF0cml4W2lkeF0gKiBsb2dnZWRNZWxCYW5kc1tuXSk7XG4gICAgfVxuICAgIG1mY2NzW2tdID0gdiAvIG51bUNvZWZmcztcbiAgfVxuXG4gIHJldHVybiBtZmNjcztcbn07IiwibW9kdWxlLmV4cG9ydHMgPSBmdW5jdGlvbihidWZmZXJTaXplLCBtKSB7XG4gIHZhciBsb3VkbmVzcyA9IG0uZmVhdHVyZUV4dHJhY3RvcnMubG91ZG5lc3MoYnVmZmVyU2l6ZSwgbSk7XG4gIHZhciBzcGVjID0gbG91ZG5lc3Muc3BlY2lmaWM7XG4gIHZhciBvdXRwdXQgPSAwO1xuXG4gIGZvciAodmFyIGkgPSAwOyBpIDwgc3BlYy5sZW5ndGg7IGkrKykge1xuICAgIGlmIChpIDwgMTUpIHtcbiAgICAgIG91dHB1dCArPSAoaSArIDEpICogc3BlY1tpICsgMV07XG4gICAgfSBlbHNlIHtcbiAgICAgIG91dHB1dCArPSAwLjA2NiAqIE1hdGguZXhwKDAuMTcxICogKGkgKyAxKSk7XG4gICAgfVxuICB9XG4gIG91dHB1dCAqPSAwLjExIC8gbG91ZG5lc3MudG90YWw7XG5cbiAgcmV0dXJuIG91dHB1dDtcbn07IiwibW9kdWxlLmV4cG9ydHMgPSBmdW5jdGlvbihidWZmZXJTaXplLCBtKSB7XG4gIHZhciBsb3VkbmVzcyA9IG0uZmVhdHVyZUV4dHJhY3RvcnMubG91ZG5lc3MoYnVmZmVyU2l6ZSwgbSk7XG5cbiAgdmFyIG1heCA9IDA7XG4gIGZvciAodmFyIGkgPSAwOyBpIDwgbG91ZG5lc3Muc3BlY2lmaWMubGVuZ3RoOyBpKyspIHtcbiAgICBpZiAobG91ZG5lc3Muc3BlY2lmaWNbaV0gPiBtYXgpIHtcbiAgICAgIG1heCA9IGxvdWRuZXNzLnNwZWNpZmljW2ldO1xuICAgIH1cbiAgfVxuXG4gIHZhciBzcHJlYWQgPSBNYXRoLnBvdygobG91ZG5lc3MudG90YWwgLSBtYXgpIC8gbG91ZG5lc3MudG90YWwsIDIpO1xuXG4gIHJldHVybiBzcHJlYWQ7XG59OyIsIm1vZHVsZS5leHBvcnRzID0gZnVuY3Rpb24oYnVmZmVyU2l6ZSwgbSkge1xuICB2YXIgcG93ZXJTcGVjdHJ1bSA9IG5ldyBGbG9hdDMyQXJyYXkobS5hbXBTcGVjdHJ1bS5sZW5ndGgpO1xuICBmb3IgKHZhciBpID0gMDsgaSA8IHBvd2VyU3BlY3RydW0ubGVuZ3RoOyBpKyspIHtcbiAgICBwb3dlclNwZWN0cnVtW2ldID0gTWF0aC5wb3cobS5hbXBTcGVjdHJ1bVtpXSwgMik7XG4gIH1cbiAgcmV0dXJuIHBvd2VyU3BlY3RydW07XG59OyIsIm1vZHVsZS5leHBvcnRzID0gZnVuY3Rpb24oYnVmZmVyU2l6ZSwgbSkge1xuXG4gIHZhciBybXMgPSAwO1xuICBmb3IgKHZhciBpID0gMDsgaSA8IG0uc2lnbmFsLmxlbmd0aDsgaSsrKSB7XG4gICAgcm1zICs9IE1hdGgucG93KG0uc2lnbmFsW2ldLCAyKTtcbiAgfVxuICBybXMgPSBybXMgLyBtLnNpZ25hbC5sZW5ndGg7XG4gIHJtcyA9IE1hdGguc3FydChybXMpO1xuXG4gIHJldHVybiBybXM7XG59O1xuIiwibW9kdWxlLmV4cG9ydHMgPSBmdW5jdGlvbihidWZmZXJTaXplLCBtKSB7XG4gIHJldHVybiDCtSgxLCBtLmFtcFNwZWN0cnVtKTtcbn07IiwibW9kdWxlLmV4cG9ydHMgPSBmdW5jdGlvbihidWZmZXJTaXplLCBtKSB7XG4gIHZhciBhbXBzcGVjID0gbS5hbXBTcGVjdHJ1bTtcbiAgdmFyIG51bWVyYXRvciA9IDA7XG4gIHZhciBkZW5vbWluYXRvciA9IDA7XG4gIGZvciAodmFyIGkgPSAwOyBpIDwgYW1wc3BlYy5sZW5ndGg7IGkrKykge1xuICAgIG51bWVyYXRvciArPSBNYXRoLmxvZyhhbXBzcGVjW2ldKTtcbiAgICBkZW5vbWluYXRvciArPSBhbXBzcGVjW2ldO1xuICB9XG4gIHJldHVybiBNYXRoLmV4cChudW1lcmF0b3IgLyBhbXBzcGVjLmxlbmd0aCkgKiBhbXBzcGVjLmxlbmd0aCAvIGRlbm9taW5hdG9yO1xufTsiLCJtb2R1bGUuZXhwb3J0cyA9IGZ1bmN0aW9uKGJ1ZmZlclNpemUsIG0pIHtcbiAgdmFyIGFtcHNwZWMgPSBtLmFtcFNwZWN0cnVtO1xuICB2YXIgwrUxID0gwrUoMSwgYW1wc3BlYyk7XG4gIHZhciDCtTIgPSDCtSgyLCBhbXBzcGVjKTtcbiAgdmFyIMK1MyA9IMK1KDMsIGFtcHNwZWMpO1xuICB2YXIgwrU0ID0gwrUoNCwgYW1wc3BlYyk7XG4gIHZhciBudW1lcmF0b3IgPSAtMyAqIE1hdGgucG93KMK1MSwgNCkgKyA2ICogwrUxICogwrUyIC0gNCAqIMK1MSAqIMK1MyArIMK1NDtcbiAgdmFyIGRlbm9taW5hdG9yID0gTWF0aC5wb3coTWF0aC5zcXJ0KMK1MiAtIE1hdGgucG93KMK1MSwgMikpLCA0KTtcbiAgcmV0dXJuIG51bWVyYXRvciAvIGRlbm9taW5hdG9yO1xufTsiLCJtb2R1bGUuZXhwb3J0cyA9IGZ1bmN0aW9uKGJ1ZmZlclNpemUsIG0pIHtcbiAgdmFyIGFtcHNwZWMgPSBtLmFtcFNwZWN0cnVtO1xuICAvL2NhbGN1bGF0ZSBueXF1aXN0IGJpblxuICB2YXIgbnlxQmluID0gbS5hdWRpb0NvbnRleHQuc2FtcGxlUmF0ZSAvICgyICogKGFtcHNwZWMubGVuZ3RoIC0gMSkpO1xuICB2YXIgZWMgPSAwO1xuICBmb3IgKHZhciBpID0gMDsgaSA8IGFtcHNwZWMubGVuZ3RoOyBpKyspIHtcbiAgICBlYyArPSBhbXBzcGVjW2ldO1xuICB9XG4gIHZhciB0aHJlc2hvbGQgPSAwLjk5ICogZWM7XG4gIHZhciBuID0gYW1wc3BlYy5sZW5ndGggLSAxO1xuICB3aGlsZSAoZWMgPiB0aHJlc2hvbGQgJiYgbiA+PSAwKSB7XG4gICAgZWMgLT0gYW1wc3BlY1tuXTtcbiAgICAtLW47XG4gIH1cbiAgcmV0dXJuIChuICsgMSkgKiBueXFCaW47XG59OyIsIm1vZHVsZS5leHBvcnRzID0gZnVuY3Rpb24oYnVmZmVyU2l6ZSwgbSwgc3BlY3RydW0pIHtcbiAgdmFyIGFtcHNwZWMgPSBtLmFtcFNwZWN0cnVtO1xuICB2YXIgwrUxID0gwrUoMSwgYW1wc3BlYyk7XG4gIHZhciDCtTIgPSDCtSgyLCBhbXBzcGVjKTtcbiAgdmFyIMK1MyA9IMK1KDMsIGFtcHNwZWMpO1xuICB2YXIgbnVtZXJhdG9yID0gMiAqIE1hdGgucG93KMK1MSwgMykgLSAzICogwrUxICogwrUyICsgwrUzO1xuICB2YXIgZGVub21pbmF0b3IgPSBNYXRoLnBvdyhNYXRoLnNxcnQowrUyIC0gTWF0aC5wb3cowrUxLCAyKSksIDMpO1xuICByZXR1cm4gbnVtZXJhdG9yIC8gZGVub21pbmF0b3I7XG59OyIsIm1vZHVsZS5leHBvcnRzID0gZnVuY3Rpb24oYnVmZmVyU2l6ZSwgbSkge1xuICAvL2xpbmVhciByZWdyZXNzaW9uXG4gIHZhciBhbXBTdW0gPSAwO1xuICB2YXIgZnJlcVN1bSA9IDA7XG4gIHZhciBmcmVxcyA9IG5ldyBGbG9hdDMyQXJyYXkobS5hbXBTcGVjdHJ1bS5sZW5ndGgpO1xuICB2YXIgcG93RnJlcVN1bSA9IDA7XG4gIHZhciBhbXBGcmVxU3VtID0gMDtcblxuICBmb3IgKHZhciBpID0gMDsgaSA8IG0uYW1wU3BlY3RydW0ubGVuZ3RoOyBpKyspIHtcbiAgICBhbXBTdW0gKz0gbS5hbXBTcGVjdHJ1bVtpXTtcbiAgICB2YXIgY3VyRnJlcSA9IGkgKiBtLmF1ZGlvQ29udGV4dC5zYW1wbGVSYXRlIC8gYnVmZmVyU2l6ZTtcbiAgICBmcmVxc1tpXSA9IGN1ckZyZXE7XG4gICAgcG93RnJlcVN1bSArPSBjdXJGcmVxICogY3VyRnJlcTtcbiAgICBmcmVxU3VtICs9IGN1ckZyZXE7XG4gICAgYW1wRnJlcVN1bSArPSBjdXJGcmVxICogbS5hbXBTcGVjdHJ1bVtpXTtcbiAgfVxuICByZXR1cm4gKG0uYW1wU3BlY3RydW0ubGVuZ3RoICogYW1wRnJlcVN1bSAtIGZyZXFTdW0gKiBhbXBTdW0pIC8gKGFtcFN1bSAqIChwb3dGcmVxU3VtIC0gTWF0aC5wb3coZnJlcVN1bSwgMikpKTtcbn07IiwibW9kdWxlLmV4cG9ydHMgPSBmdW5jdGlvbihidWZmZXJTaXplLCBtKSB7XG4gIHZhciBhbXBzcGVjID0gbS5hbXBTcGVjdHJ1bTtcbiAgcmV0dXJuIE1hdGguc3FydCjCtSgyLCBhbXBzcGVjKSAtIE1hdGgucG93KMK1KDEsIGFtcHNwZWMpLCAyKSk7XG59OyIsIm1vZHVsZS5leHBvcnRzID0gZnVuY3Rpb24oYnVmZmVyU2l6ZSwgbSkge1xuICB2YXIgemNyID0gMDtcbiAgZm9yICh2YXIgaSA9IDA7IGkgPCBtLnNpZ25hbC5sZW5ndGg7IGkrKykge1xuICAgIGlmICgobS5zaWduYWxbaV0gPj0gMCAmJiBtLnNpZ25hbFtpICsgMV0gPCAwKSB8fCAobS5zaWduYWxbaV0gPCAwICYmIG0uc2lnbmFsW2kgKyAxXSA+PSAwKSkge1xuICAgICAgemNyKys7XG4gICAgfVxuICB9XG4gIHJldHVybiB6Y3I7XG59OyIsIi8vIE1leWRhIEphdmFzY3JpcHQgRFNQIGxpYnJhcnlcblxudmFyIENvbXBsZXhBcnJheSA9IHJlcXVpcmUoJy4uL2xpYi9qc2ZmdC9jb21wbGV4X2FycmF5JykuQ29tcGxleEFycmF5O1xuXG4vLyBtb2RpZmllcyBDb21wbGV4QXJyYXlcbnZhciBmZnQgPSByZXF1aXJlKCcuLi9saWIvanNmZnQvZmZ0Jyk7XG5cbm1vZHVsZS5leHBvcnRzID0gZnVuY3Rpb24oYXVkaW9Db250ZXh0LHNyYyxidWZTaXplLGNhbGxiYWNrKXtcblx0XG5cdC8vSSBhbSBteXNlbGZcblx0dmFyIHNlbGYgPSB0aGlzO1xuXHRcblx0c2VsZi5mZWF0dXJlRXh0cmFjdG9ycyA9IHJlcXVpcmUoJy4vZXh0cmFjdG9ycycpO1xuXG5cdC8vZGVmYXVsdCBidWZmZXIgc2l6ZVxuXHR2YXIgYnVmZmVyU2l6ZSA9IGJ1ZlNpemUgPyBidWZTaXplIDogMjU2O1xuXG5cdC8vaW5pdGlhbCBzb3VyY2Vcblx0dmFyIHNvdXJjZSA9IHNyYztcblxuXHQvL2NhbGxiYWNrIGNvbnRyb2xsZXJzXG5cdHZhciBFWFRSQUNUSU9OX1NUQVJURUQgPSBmYWxzZTtcblx0dmFyIF9mZWF0dXJlc1RvRXh0cmFjdDtcblxuXHQvL3V0aWxpdGllc1xuXHR2YXIgwrUgPSBmdW5jdGlvbihpLCBhbXBsaXR1ZGVTcGVjdCl7XG5cdFx0dmFyIG51bWVyYXRvciA9IDA7XG5cdFx0dmFyIGRlbm9taW5hdG9yID0gMDtcblx0XHRmb3IodmFyIGsgPSAwOyBrIDwgYW1wbGl0dWRlU3BlY3QubGVuZ3RoOyBrKyspe1xuXHRcdFx0bnVtZXJhdG9yICs9IE1hdGgucG93KGssaSkqTWF0aC5hYnMoYW1wbGl0dWRlU3BlY3Rba10pO1xuXHRcdFx0ZGVub21pbmF0b3IgKz0gYW1wbGl0dWRlU3BlY3Rba107XG5cdFx0fVxuXHRcdHJldHVybiBudW1lcmF0b3IvZGVub21pbmF0b3I7XG5cdH07XG5cblx0dmFyIGlzUG93ZXJPZlR3byA9IGZ1bmN0aW9uKG51bSkge1xuXHRcdHdoaWxlICgoKG51bSAlIDIpID09IDApICYmIG51bSA+IDEpIHtcblx0XHRcdG51bSAvPSAyO1xuXHRcdH1cblx0XHRyZXR1cm4gKG51bSA9PSAxKTtcblx0fTtcblxuXHQvL2luaXRpbGl6ZSBiYXJrIHNjYWxlICh0byBiZSB1c2VkIGluIG1vc3QgcGVyY2VwdHVhbCBmZWF0dXJlcykuXG5cdHNlbGYuYmFya1NjYWxlID0gbmV3IEZsb2F0MzJBcnJheShidWZTaXplKTtcblxuXHRmb3IodmFyIGkgPSAwOyBpIDwgc2VsZi5iYXJrU2NhbGUubGVuZ3RoOyBpKyspe1xuXHRcdHNlbGYuYmFya1NjYWxlW2ldID0gaSphdWRpb0NvbnRleHQuc2FtcGxlUmF0ZS8oYnVmU2l6ZSk7XG5cdFx0c2VsZi5iYXJrU2NhbGVbaV0gPSAxMypNYXRoLmF0YW4oc2VsZi5iYXJrU2NhbGVbaV0vMTMxNS44KSArIDMuNSogTWF0aC5hdGFuKE1hdGgucG93KChzZWxmLmJhcmtTY2FsZVtpXS83NTE4KSwyKSk7XG5cdH1cblxuXHQvL1dJTkRPV0lOR1xuXHQvL3NldCBkZWZhdWx0XG5cdHNlbGYud2luZG93aW5nRnVuY3Rpb24gPSBcImhhbm5pbmdcIjtcblxuXHQvL2NyZWF0ZSB3aW5kb3dzXG5cdHNlbGYuaGFubmluZyA9IG5ldyBGbG9hdDMyQXJyYXkoYnVmU2l6ZSk7XG5cdGZvciAodmFyIGkgPSAwOyBpIDwgYnVmU2l6ZTsgaSsrKSB7XG5cdFx0Ly9BY2NvcmRpbmcgdG8gdGhlIFIgZG9jdW1lbnRhdGlvbiBodHRwOi8vcmdtLm9nYWxhYi5uZXQvUkdNL1JfcmRmaWxlP2Y9R0VORUFyZWFkL21hbi9oYW5uaW5nLndpbmRvdy5SZCZkPVJfQ0Ncblx0XHRzZWxmLmhhbm5pbmdbaV0gPSAwLjUgLSAwLjUqTWF0aC5jb3MoMipNYXRoLlBJKmkvKGJ1ZlNpemUtMSkpO1xuXHR9XG5cblx0c2VsZi5oYW1taW5nID0gbmV3IEZsb2F0MzJBcnJheShidWZTaXplKTtcblx0Zm9yICh2YXIgaSA9IDA7IGkgPCBidWZTaXplOyBpKyspIHtcblx0XHQvL0FjY29yZGluZyB0byBodHRwOi8vdWsubWF0aHdvcmtzLmNvbS9oZWxwL3NpZ25hbC9yZWYvaGFtbWluZy5odG1sXG5cdFx0c2VsZi5oYW1taW5nW2ldID0gMC41NCAtIDAuNDYqTWF0aC5jb3MoMipNYXRoLlBJKihpL2J1ZlNpemUtMSkpO1xuXHR9XG5cblx0Ly9VTkZJTklTSEVEIC0gYmxhY2ttYW4gd2luZG93IGltcGxlbWVudGF0aW9uXG5cblx0LypzZWxmLmJsYWNrbWFuID0gbmV3IEZsb2F0MzJBcnJheShidWZTaXplKTtcblx0Ly9BY2NvcmRpbmcgdG8gaHR0cDovL3VrLm1hdGh3b3Jrcy5jb20vaGVscC9zaWduYWwvcmVmL2JsYWNrbWFuLmh0bWxcblx0Ly9maXJzdCBoYWxmIG9mIHRoZSB3aW5kb3dcblx0Zm9yICh2YXIgaSA9IDA7IGkgPCAoYnVmU2l6ZSAlIDIpID8gKGJ1ZlNpemUrMSkvMiA6IGJ1ZlNpemUvMjsgaSsrKSB7XG5cdFx0c2VsZi5ibGFja21hbltpXSA9IDAuNDIgLSAwLjUqTWF0aC5jb3MoMipNYXRoLlBJKmkvKGJ1ZlNpemUtMSkpICsgMC4wOCpNYXRoLmNvcyg0Kk1hdGguUEkqaS8oYnVmU2l6ZS0xKSk7XG5cdH1cblx0Ly9zZWNvbmQgaGFsZiBvZiB0aGUgd2luZG93XG5cdGZvciAodmFyIGkgPSBidWZTaXplLzI7IGkgPiAwOyBpLS0pIHtcblx0XHRzZWxmLmJsYWNrbWFuW2J1ZlNpemUgLSBpXSA9IHNlbGYuYmxhY2ttYW5baV07XG5cdH0qL1xuXG5cdHNlbGYud2luZG93aW5nID0gZnVuY3Rpb24oc2lnLCB0eXBlKXtcblx0XHR2YXIgd2luZG93ZWQgPSBuZXcgRmxvYXQzMkFycmF5KHNpZy5sZW5ndGgpO1xuXHRcdHZhciBpICxsZW4gPSBzaWcubGVuZ3RoO1xuXG5cdFx0aWYgKHR5cGUgPT0gXCJoYW5uaW5nXCIpIHtcblx0XHRcdGZvciAoaSA9IDA7IGkgPCBsZW47IGkrKykge1xuXHRcdFx0XHR3aW5kb3dlZFtpXSA9IHNpZ1tpXSpzZWxmLmhhbm5pbmdbaV07XG5cdFx0XHR9XG5cdFx0fVxuXHRcdGVsc2UgaWYgKHR5cGUgPT0gXCJoYW1taW5nXCIpIHtcblx0XHRcdGZvciAoaSA9IDA7IGkgPCBsZW47IGkrKykge1xuXHRcdFx0XHR3aW5kb3dlZFtpXSA9IHNpZ1tpXSpzZWxmLmhhbW1pbmdbaV07XG5cdFx0XHR9XG5cdFx0fVxuXHRcdGVsc2UgaWYgKHR5cGUgPT0gXCJibGFja21hblwiKSB7XG5cdFx0XHRmb3IgKGkgPSAwOyBpIDwgbGVuOyBpKyspIHtcblx0XHRcdFx0d2luZG93ZWRbaV0gPSBzaWdbaV0qc2VsZi5ibGFja21hbltpXTtcblx0XHRcdH1cblx0XHR9XG5cblx0XHRyZXR1cm4gd2luZG93ZWQ7XG5cdH07XG5cblx0Ly9zb3VyY2Ugc2V0dGVyIG1ldGhvZFxuXHRzZWxmLnNldFNvdXJjZSA9IGZ1bmN0aW9uKF9zcmMpIHtcblx0XHRzb3VyY2UgPSBfc3JjO1xuXHRcdHNvdXJjZS5jb25uZWN0KHdpbmRvdy5zcG4pO1xuXHR9O1xuXG5cblxuXHRpZiAoaXNQb3dlck9mVHdvKGJ1ZmZlclNpemUpICYmIGF1ZGlvQ29udGV4dCkge1xuXHRcdFx0c2VsZi5mZWF0dXJlSW5mbyA9IHtcblx0XHRcdFx0XCJidWZmZXJcIjoge1xuXHRcdFx0XHRcdFwidHlwZVwiOiBcImFycmF5XCJcblx0XHRcdFx0fSxcblx0XHRcdFx0XCJybXNcIjoge1xuXHRcdFx0XHRcdFwidHlwZVwiOiBcIm51bWJlclwiXG5cdFx0XHRcdH0sXG5cdFx0XHRcdFwiZW5lcmd5XCI6IHtcblx0XHRcdFx0XHRcInR5cGVcIjogXCJudW1iZXJcIlxuXHRcdFx0XHR9LFxuXHRcdFx0XHRcInpjclwiOiB7XG5cdFx0XHRcdFx0XCJ0eXBlXCI6IFwibnVtYmVyXCJcblx0XHRcdFx0fSxcblx0XHRcdFx0XCJjb21wbGV4U3BlY3RydW1cIjoge1xuXHRcdFx0XHRcdFwidHlwZVwiOiBcIm11bHRpcGxlQXJyYXlzXCIsXG5cdFx0XHRcdFx0XCJhcnJheU5hbWVzXCI6IHtcblx0XHRcdFx0XHRcdFwiMVwiOiBcInJlYWxcIixcblx0XHRcdFx0XHRcdFwiMlwiOiBcImltYWdcIlxuXHRcdFx0XHRcdH1cblx0XHRcdFx0fSxcblx0XHRcdFx0XCJhbXBsaXR1ZGVTcGVjdHJ1bVwiOiB7XG5cdFx0XHRcdFx0XCJ0eXBlXCI6IFwiYXJyYXlcIlxuXHRcdFx0XHR9LFxuXHRcdFx0XHRcInBvd2VyU3BlY3RydW1cIjoge1xuXHRcdFx0XHRcdFwidHlwZVwiOiBcImFycmF5XCJcblx0XHRcdFx0fSxcblx0XHRcdFx0XCJzcGVjdHJhbENlbnRyb2lkXCI6IHtcblx0XHRcdFx0XHRcInR5cGVcIjogXCJudW1iZXJcIlxuXHRcdFx0XHR9LFxuXHRcdFx0XHRcInNwZWN0cmFsRmxhdG5lc3NcIjoge1xuXHRcdFx0XHRcdFwidHlwZVwiOiBcIm51bWJlclwiXG5cdFx0XHRcdH0sXG5cdFx0XHRcdFwic3BlY3RyYWxTbG9wZVwiOiB7XG5cdFx0XHRcdFx0XCJ0eXBlXCI6IFwibnVtYmVyXCJcblx0XHRcdFx0fSxcblx0XHRcdFx0XCJzcGVjdHJhbFJvbGxvZmZcIjoge1xuXHRcdFx0XHRcdFwidHlwZVwiOiBcIm51bWJlclwiXG5cdFx0XHRcdH0sXG5cdFx0XHRcdFwic3BlY3RyYWxTcHJlYWRcIjoge1xuXHRcdFx0XHRcdFwidHlwZVwiOiBcIm51bWJlclwiXG5cdFx0XHRcdH0sXG5cdFx0XHRcdFwic3BlY3RyYWxTa2V3bmVzc1wiOiB7XG5cdFx0XHRcdFx0XCJ0eXBlXCI6IFwibnVtYmVyXCJcblx0XHRcdFx0fSxcblx0XHRcdFx0XCJzcGVjdHJhbEt1cnRvc2lzXCI6IHtcblx0XHRcdFx0XHRcInR5cGVcIjogXCJudW1iZXJcIlxuXHRcdFx0XHR9LFxuXHRcdFx0XHRcImxvdWRuZXNzXCI6IHtcblx0XHRcdFx0XHRcInR5cGVcIjogXCJtdWx0aXBsZUFycmF5c1wiLFxuXHRcdFx0XHRcdFwiYXJyYXlOYW1lc1wiOiB7XG5cdFx0XHRcdFx0XHRcIjFcIjogXCJ0b3RhbFwiLFxuXHRcdFx0XHRcdFx0XCIyXCI6IFwic3BlY2lmaWNcIlxuXHRcdFx0XHRcdH1cblx0XHRcdFx0fSxcblx0XHRcdFx0XCJwZXJjZXB0dWFsU3ByZWFkXCI6IHtcblx0XHRcdFx0XHRcInR5cGVcIjogXCJudW1iZXJcIlxuXHRcdFx0XHR9LFxuXHRcdFx0XHRcInBlcmNlcHR1YWxTaGFycG5lc3NcIjoge1xuXHRcdFx0XHRcdFwidHlwZVwiOiBcIm51bWJlclwiXG5cdFx0XHRcdH0sXG5cdFx0XHRcdFwibWZjY1wiOiB7XG5cdFx0XHRcdFx0XCJ0eXBlXCI6IFwiYXJyYXlcIlxuXHRcdFx0XHR9XG5cdFx0XHR9O1xuXG5cdFx0XHQvL2NyZWF0ZSBjb21wbGV4YXJyYXkgdG8gaG9sZCB0aGUgc3BlY3RydW1cblx0XHRcdHZhciBkYXRhID0gbmV3IENvbXBsZXhBcnJheShidWZmZXJTaXplKTtcblx0XHRcdFxuXHRcdFx0Ly90cmFuc2Zvcm1cblx0XHRcdHZhciBzcGVjID0gZGF0YS5GRlQoKTtcblx0XHRcdC8vYXNzaWduIHRvIG1leWRhXG5cdFx0XHRzZWxmLmNvbXBsZXhTcGVjdHJ1bSA9IHNwZWM7XG5cdFx0XHRzZWxmLmFtcFNwZWN0cnVtID0gbmV3IEZsb2F0MzJBcnJheShidWZmZXJTaXplLzIpO1xuXG5cdFx0XHQvL2NyZWF0ZSBub2Rlc1xuXHRcdFx0d2luZG93LnNwbiA9IGF1ZGlvQ29udGV4dC5jcmVhdGVTY3JpcHRQcm9jZXNzb3IoYnVmZmVyU2l6ZSwxLDEpO1xuXHRcdFx0c3BuLmNvbm5lY3QoYXVkaW9Db250ZXh0LmRlc3RpbmF0aW9uKTtcblxuXHRcdFx0d2luZG93LnNwbi5vbmF1ZGlvcHJvY2VzcyA9IGZ1bmN0aW9uKGUpIHtcblx0XHRcdFx0Ly90aGlzIGlzIHRvIG9idGFpbiB0aGUgY3VycmVudCBhbXBsaXR1ZGUgc3BlY3RydW1cblx0XHRcdFx0dmFyIGlucHV0RGF0YSA9IGUuaW5wdXRCdWZmZXIuZ2V0Q2hhbm5lbERhdGEoMCk7XG5cdFx0XHRcdHNlbGYuc2lnbmFsID0gaW5wdXREYXRhO1xuXHRcdFx0XHR2YXIgd2luZG93ZWRTaWduYWwgPSBzZWxmLndpbmRvd2luZyhzZWxmLnNpZ25hbCwgc2VsZi53aW5kb3dpbmdGdW5jdGlvbik7XG5cblx0XHRcdFx0Ly9tYXAgdGltZSBkb21haW5cblx0XHRcdFx0ZGF0YS5tYXAoZnVuY3Rpb24odmFsdWUsIGksIG4pIHtcblx0XHRcdFx0XHR2YWx1ZS5yZWFsID0gd2luZG93ZWRTaWduYWxbaV07XG5cdFx0XHRcdH0pO1xuXG5cdFx0XHRcdC8vY2FsY3VsYXRlIGFtcGxpdHVkZVxuXHRcdFx0XHRmb3IgKHZhciBpID0gMDsgaSA8IGJ1ZmZlclNpemUvMjsgaSsrKSB7XG5cdFx0XHRcdFx0c2VsZi5hbXBTcGVjdHJ1bVtpXSA9IE1hdGguc3FydChNYXRoLnBvdyhzcGVjLnJlYWxbaV0sMikgKyBNYXRoLnBvdyhzcGVjLmltYWdbaV0sMikpO1xuXHRcdFx0XHR9XG5cblx0XHRcdFx0Ly9jYWxsIGNhbGxiYWNrIGlmIGFwcGxpY2FibGVcblx0XHRcdFx0aWYgKHR5cGVvZiBjYWxsYmFjayA9PT0gXCJmdW5jdGlvblwiICYmIEVYVFJBQ1RJT05fU1RBUlRFRCkge1xuXHRcdFx0XHRcdGNhbGxiYWNrKHNlbGYuZ2V0KF9mZWF0dXJlc1RvRXh0cmFjdCkpO1xuXHRcdFx0XHR9XG5cblx0XHRcdH07XG5cblx0XHRcdHNlbGYuc3RhcnQgPSBmdW5jdGlvbihmZWF0dXJlcykge1xuXHRcdFx0XHRfZmVhdHVyZXNUb0V4dHJhY3QgPSBmZWF0dXJlcztcblx0XHRcdFx0RVhUUkFDVElPTl9TVEFSVEVEID0gdHJ1ZTtcblx0XHRcdH07XG5cblx0XHRcdHNlbGYuc3RvcCA9IGZ1bmN0aW9uKCkge1xuXHRcdFx0XHRFWFRSQUNUSU9OX1NUQVJURUQgPSBmYWxzZTtcblx0XHRcdH07XG5cblx0XHRcdHNlbGYuYXVkaW9Db250ZXh0ID0gYXVkaW9Db250ZXh0O1xuXG5cdFx0XHRzZWxmLmdldCA9IGZ1bmN0aW9uKGZlYXR1cmUpIHtcblx0XHRcdFx0aWYodHlwZW9mIGZlYXR1cmUgPT09IFwib2JqZWN0XCIpe1xuXHRcdFx0XHRcdHZhciByZXN1bHRzID0ge307XG5cdFx0XHRcdFx0Zm9yICh2YXIgeCA9IDA7IHggPCBmZWF0dXJlLmxlbmd0aDsgeCsrKXtcblx0XHRcdFx0XHRcdHRyeXtcblx0XHRcdFx0XHRcdFx0cmVzdWx0c1tmZWF0dXJlW3hdXSA9IChzZWxmLmZlYXR1cmVFeHRyYWN0b3JzW2ZlYXR1cmVbeF1dKGJ1ZmZlclNpemUsIHNlbGYpKTtcblx0XHRcdFx0XHRcdH0gY2F0Y2ggKGUpe1xuXHRcdFx0XHRcdFx0XHRjb25zb2xlLmVycm9yKGUpO1xuXHRcdFx0XHRcdFx0fVxuXHRcdFx0XHRcdH1cblx0XHRcdFx0XHRyZXR1cm4gcmVzdWx0cztcblx0XHRcdFx0fVxuXHRcdFx0XHRlbHNlIGlmICh0eXBlb2YgZmVhdHVyZSA9PT0gXCJzdHJpbmdcIil7XG5cdFx0XHRcdFx0cmV0dXJuIHNlbGYuZmVhdHVyZUV4dHJhY3RvcnNbZmVhdHVyZV0oYnVmZmVyU2l6ZSwgc2VsZik7XG5cdFx0XHRcdH1cblx0XHRcdFx0ZWxzZXtcblx0XHRcdFx0XHR0aHJvdyBcIkludmFsaWQgRmVhdHVyZSBGb3JtYXRcIjtcblx0XHRcdFx0fVxuXHRcdFx0fTtcblx0XHRcdHNvdXJjZS5jb25uZWN0KHdpbmRvdy5zcG4sIDAsIDApO1xuXHRcdFx0cmV0dXJuIHNlbGY7XG5cdH1cblx0ZWxzZSB7XG5cdFx0Ly9oYW5kbGUgZXJyb3JzXG5cdFx0aWYgKHR5cGVvZiBhdWRpb0NvbnRleHQgPT0gXCJ1bmRlZmluZWRcIikge1xuXHRcdFx0dGhyb3cgXCJBdWRpb0NvbnRleHQgd2Fzbid0IHNwZWNpZmllZDogTWV5ZGEgd2lsbCBub3QgcnVuLlwiO1xuXHRcdH1cblx0XHRlbHNlIHtcblx0XHRcdHRocm93IFwiQnVmZmVyIHNpemUgaXMgbm90IGEgcG93ZXIgb2YgdHdvOiBNZXlkYSB3aWxsIG5vdCBydW4uXCI7XG5cdFx0fVxuXHR9XG59OyJdfQ==

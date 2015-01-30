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
},{"../lib/jsfft/complex_array":1,"../lib/jsfft/fft":2,"./extractors":7}]},{},[22])(22)
});
//# sourceMappingURL=data:application/json;base64,eyJ2ZXJzaW9uIjozLCJzb3VyY2VzIjpbIi4uLy4uLy4uLy4uL3Vzci9sb2NhbC9saWIvbm9kZV9tb2R1bGVzL2Jyb3dzZXJpZnkvbm9kZV9tb2R1bGVzL2Jyb3dzZXItcGFjay9fcHJlbHVkZS5qcyIsImxpYi9qc2ZmdC9jb21wbGV4X2FycmF5LmpzIiwibGliL2pzZmZ0L2ZmdC5qcyIsInNyYy9leHRyYWN0b3JzL2FtcGxpdHVkZVNwZWN0cnVtLmpzIiwic3JjL2V4dHJhY3RvcnMvYnVmZmVyLmpzIiwic3JjL2V4dHJhY3RvcnMvY29tcGxleFNwZWN0cnVtLmpzIiwic3JjL2V4dHJhY3RvcnMvZW5lcmd5LmpzIiwic3JjL2V4dHJhY3RvcnMvaW5kZXguanMiLCJzcmMvZXh0cmFjdG9ycy9sb3VkbmVzcy5qcyIsInNyYy9leHRyYWN0b3JzL21mY2MuanMiLCJzcmMvZXh0cmFjdG9ycy9wZXJjZXB0dWFsU2hhcnBuZXNzLmpzIiwic3JjL2V4dHJhY3RvcnMvcGVyY2VwdHVhbFNwcmVhZC5qcyIsInNyYy9leHRyYWN0b3JzL3Bvd2VyU3BlY3RydW0uanMiLCJzcmMvZXh0cmFjdG9ycy9ybXMuanMiLCJzcmMvZXh0cmFjdG9ycy9zcGVjdHJhbENlbnRyb2lkLmpzIiwic3JjL2V4dHJhY3RvcnMvc3BlY3RyYWxGbGF0bmVzcy5qcyIsInNyYy9leHRyYWN0b3JzL3NwZWN0cmFsS3VydG9zaXMuanMiLCJzcmMvZXh0cmFjdG9ycy9zcGVjdHJhbFJvbGxvZmYuanMiLCJzcmMvZXh0cmFjdG9ycy9zcGVjdHJhbFNrZXduZXNzLmpzIiwic3JjL2V4dHJhY3RvcnMvc3BlY3RyYWxTbG9wZS5qcyIsInNyYy9leHRyYWN0b3JzL3NwZWN0cmFsU3ByZWFkLmpzIiwic3JjL2V4dHJhY3RvcnMvemNyLmpzIiwic3JjL21leWRhLmpzIl0sIm5hbWVzIjpbXSwibWFwcGluZ3MiOiJBQUFBO0FDQUE7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBOztBQ3BIQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTs7QUNqT0E7QUFDQTtBQUNBOztBQ0ZBO0FBQ0E7QUFDQTs7QUNGQTtBQUNBO0FBQ0E7O0FDRkE7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7O0FDTkE7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTs7QUNuQkE7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7O0FDN0NBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTs7QUMvRkE7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7O0FDZkE7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTs7QUNiQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTs7QUNOQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7O0FDWEE7QUFDQTtBQUNBOztBQ0ZBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBOztBQ1RBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBOztBQ1RBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBOztBQ2ZBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTs7QUNSQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7O0FDakJBO0FBQ0E7QUFDQTtBQUNBOztBQ0hBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTs7QUNSQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBIiwiZmlsZSI6ImdlbmVyYXRlZC5qcyIsInNvdXJjZVJvb3QiOiIiLCJzb3VyY2VzQ29udGVudCI6WyIoZnVuY3Rpb24gZSh0LG4scil7ZnVuY3Rpb24gcyhvLHUpe2lmKCFuW29dKXtpZighdFtvXSl7dmFyIGE9dHlwZW9mIHJlcXVpcmU9PVwiZnVuY3Rpb25cIiYmcmVxdWlyZTtpZighdSYmYSlyZXR1cm4gYShvLCEwKTtpZihpKXJldHVybiBpKG8sITApO3ZhciBmPW5ldyBFcnJvcihcIkNhbm5vdCBmaW5kIG1vZHVsZSAnXCIrbytcIidcIik7dGhyb3cgZi5jb2RlPVwiTU9EVUxFX05PVF9GT1VORFwiLGZ9dmFyIGw9bltvXT17ZXhwb3J0czp7fX07dFtvXVswXS5jYWxsKGwuZXhwb3J0cyxmdW5jdGlvbihlKXt2YXIgbj10W29dWzFdW2VdO3JldHVybiBzKG4/bjplKX0sbCxsLmV4cG9ydHMsZSx0LG4scil9cmV0dXJuIG5bb10uZXhwb3J0c312YXIgaT10eXBlb2YgcmVxdWlyZT09XCJmdW5jdGlvblwiJiZyZXF1aXJlO2Zvcih2YXIgbz0wO288ci5sZW5ndGg7bysrKXMocltvXSk7cmV0dXJuIHN9KSIsIid1c2Ugc3RyaWN0JztcblxuIWZ1bmN0aW9uKGV4cG9ydHMsIHVuZGVmaW5lZCkge1xuXG4gIHZhclxuICAgIC8vIElmIHRoZSB0eXBlZCBhcnJheSBpcyB1bnNwZWNpZmllZCwgdXNlIHRoaXMuXG4gICAgRGVmYXVsdEFycmF5VHlwZSA9IEZsb2F0MzJBcnJheSxcbiAgICAvLyBTaW1wbGUgbWF0aCBmdW5jdGlvbnMgd2UgbmVlZC5cbiAgICBzcXJ0ID0gTWF0aC5zcXJ0LFxuICAgIHNxciA9IGZ1bmN0aW9uKG51bWJlcikge3JldHVybiBNYXRoLnBvdyhudW1iZXIsIDIpfSxcbiAgICAvLyBJbnRlcm5hbCBjb252ZW5pZW5jZSBjb3BpZXMgb2YgdGhlIGV4cG9ydGVkIGZ1bmN0aW9uc1xuICAgIGlzQ29tcGxleEFycmF5LFxuICAgIENvbXBsZXhBcnJheVxuXG4gIGV4cG9ydHMuaXNDb21wbGV4QXJyYXkgPSBpc0NvbXBsZXhBcnJheSA9IGZ1bmN0aW9uKG9iaikge1xuICAgIHJldHVybiBvYmogIT09IHVuZGVmaW5lZCAmJlxuICAgICAgb2JqLmhhc093blByb3BlcnR5ICE9PSB1bmRlZmluZWQgJiZcbiAgICAgIG9iai5oYXNPd25Qcm9wZXJ0eSgncmVhbCcpICYmXG4gICAgICBvYmouaGFzT3duUHJvcGVydHkoJ2ltYWcnKVxuICB9XG5cbiAgZXhwb3J0cy5Db21wbGV4QXJyYXkgPSBDb21wbGV4QXJyYXkgPSBmdW5jdGlvbihvdGhlciwgb3B0X2FycmF5X3R5cGUpe1xuICAgIGlmIChpc0NvbXBsZXhBcnJheShvdGhlcikpIHtcbiAgICAgIC8vIENvcHkgY29uc3R1Y3Rvci5cbiAgICAgIHRoaXMuQXJyYXlUeXBlID0gb3RoZXIuQXJyYXlUeXBlXG4gICAgICB0aGlzLnJlYWwgPSBuZXcgdGhpcy5BcnJheVR5cGUob3RoZXIucmVhbClcbiAgICAgIHRoaXMuaW1hZyA9IG5ldyB0aGlzLkFycmF5VHlwZShvdGhlci5pbWFnKVxuICAgIH0gZWxzZSB7XG4gICAgICB0aGlzLkFycmF5VHlwZSA9IG9wdF9hcnJheV90eXBlIHx8IERlZmF1bHRBcnJheVR5cGVcbiAgICAgIC8vIG90aGVyIGNhbiBiZSBlaXRoZXIgYW4gYXJyYXkgb3IgYSBudW1iZXIuXG4gICAgICB0aGlzLnJlYWwgPSBuZXcgdGhpcy5BcnJheVR5cGUob3RoZXIpXG4gICAgICB0aGlzLmltYWcgPSBuZXcgdGhpcy5BcnJheVR5cGUodGhpcy5yZWFsLmxlbmd0aClcbiAgICB9XG5cbiAgICB0aGlzLmxlbmd0aCA9IHRoaXMucmVhbC5sZW5ndGhcbiAgfVxuXG4gIENvbXBsZXhBcnJheS5wcm90b3R5cGUudG9TdHJpbmcgPSBmdW5jdGlvbigpIHtcbiAgICB2YXIgY29tcG9uZW50cyA9IFtdXG5cbiAgICB0aGlzLmZvckVhY2goZnVuY3Rpb24oY192YWx1ZSwgaSkge1xuICAgICAgY29tcG9uZW50cy5wdXNoKFxuICAgICAgICAnKCcgK1xuICAgICAgICBjX3ZhbHVlLnJlYWwudG9GaXhlZCgyKSArICcsJyArXG4gICAgICAgIGNfdmFsdWUuaW1hZy50b0ZpeGVkKDIpICtcbiAgICAgICAgJyknXG4gICAgICApXG4gICAgfSlcblxuICAgIHJldHVybiAnWycgKyBjb21wb25lbnRzLmpvaW4oJywnKSArICddJ1xuICB9XG5cbiAgLy8gSW4tcGxhY2UgbWFwcGVyLlxuICBDb21wbGV4QXJyYXkucHJvdG90eXBlLm1hcCA9IGZ1bmN0aW9uKG1hcHBlcikge1xuICAgIHZhclxuICAgICAgaSxcbiAgICAgIG4gPSB0aGlzLmxlbmd0aCxcbiAgICAgIC8vIEZvciBHQyBlZmZpY2llbmN5LCBwYXNzIGEgc2luZ2xlIGNfdmFsdWUgb2JqZWN0IHRvIHRoZSBtYXBwZXIuXG4gICAgICBjX3ZhbHVlID0ge31cblxuICAgIGZvciAoaSA9IDA7IGkgPCBuOyBpKyspIHtcbiAgICAgIGNfdmFsdWUucmVhbCA9IHRoaXMucmVhbFtpXVxuICAgICAgY192YWx1ZS5pbWFnID0gdGhpcy5pbWFnW2ldXG4gICAgICBtYXBwZXIoY192YWx1ZSwgaSwgbilcbiAgICAgIHRoaXMucmVhbFtpXSA9IGNfdmFsdWUucmVhbFxuICAgICAgdGhpcy5pbWFnW2ldID0gY192YWx1ZS5pbWFnXG4gICAgfVxuXG4gICAgcmV0dXJuIHRoaXNcbiAgfVxuXG4gIENvbXBsZXhBcnJheS5wcm90b3R5cGUuZm9yRWFjaCA9IGZ1bmN0aW9uKGl0ZXJhdG9yKSB7XG4gICAgdmFyXG4gICAgICBpLFxuICAgICAgbiA9IHRoaXMubGVuZ3RoLFxuICAgICAgLy8gRm9yIGNvbnNpc3RlbmN5IHdpdGggLm1hcC5cbiAgICAgIGNfdmFsdWUgPSB7fVxuXG4gICAgZm9yIChpID0gMDsgaSA8IG47IGkrKykge1xuICAgICAgY192YWx1ZS5yZWFsID0gdGhpcy5yZWFsW2ldXG4gICAgICBjX3ZhbHVlLmltYWcgPSB0aGlzLmltYWdbaV1cbiAgICAgIGl0ZXJhdG9yKGNfdmFsdWUsIGksIG4pXG4gICAgfVxuICB9XG5cbiAgQ29tcGxleEFycmF5LnByb3RvdHlwZS5jb25qdWdhdGUgPSBmdW5jdGlvbigpIHtcbiAgICByZXR1cm4gKG5ldyBDb21wbGV4QXJyYXkodGhpcykpLm1hcChmdW5jdGlvbih2YWx1ZSkge1xuICAgICAgdmFsdWUuaW1hZyAqPSAtMVxuICAgIH0pXG4gIH1cblxuICAvLyBIZWxwZXIgc28gd2UgY2FuIG1ha2UgQXJyYXlUeXBlIG9iamVjdHMgcmV0dXJuZWQgaGF2ZSBzaW1pbGFyIGludGVyZmFjZXNcbiAgLy8gICB0byBDb21wbGV4QXJyYXlzLlxuICBmdW5jdGlvbiBpdGVyYWJsZShvYmopIHtcbiAgICBpZiAoIW9iai5mb3JFYWNoKVxuICAgICAgb2JqLmZvckVhY2ggPSBmdW5jdGlvbihpdGVyYXRvcikge1xuICAgICAgICB2YXIgaSwgbiA9IHRoaXMubGVuZ3RoXG5cbiAgICAgICAgZm9yIChpID0gMDsgaSA8IG47IGkrKylcbiAgICAgICAgICBpdGVyYXRvcih0aGlzW2ldLCBpLCBuKVxuICAgICAgfVxuXG4gICAgcmV0dXJuIG9ialxuICB9XG5cbiAgQ29tcGxleEFycmF5LnByb3RvdHlwZS5tYWduaXR1ZGUgPSBmdW5jdGlvbigpIHtcbiAgICB2YXIgbWFncyA9IG5ldyB0aGlzLkFycmF5VHlwZSh0aGlzLmxlbmd0aClcblxuICAgIHRoaXMuZm9yRWFjaChmdW5jdGlvbih2YWx1ZSwgaSkge1xuICAgICAgbWFnc1tpXSA9IHNxcnQoc3FyKHZhbHVlLnJlYWwpICsgc3FyKHZhbHVlLmltYWcpKVxuICAgIH0pXG5cbiAgICAvLyBBcnJheVR5cGUgd2lsbCBub3QgbmVjZXNzYXJpbHkgYmUgaXRlcmFibGU6IG1ha2UgaXQgc28uXG4gICAgcmV0dXJuIGl0ZXJhYmxlKG1hZ3MpXG4gIH1cbn0odHlwZW9mIGV4cG9ydHMgPT09ICd1bmRlZmluZWQnICYmICh0aGlzLmNvbXBsZXhfYXJyYXkgPSB7fSkgfHwgZXhwb3J0cylcbiIsIid1c2Ugc3RyaWN0JztcblxuIWZ1bmN0aW9uKGV4cG9ydHMsIGNvbXBsZXhfYXJyYXkpIHtcblxuICB2YXJcbiAgICBDb21wbGV4QXJyYXkgPSBjb21wbGV4X2FycmF5LkNvbXBsZXhBcnJheSxcbiAgICAvLyBNYXRoIGNvbnN0YW50cyBhbmQgZnVuY3Rpb25zIHdlIG5lZWQuXG4gICAgUEkgPSBNYXRoLlBJLFxuICAgIFNRUlQxXzIgPSBNYXRoLlNRUlQxXzIsXG4gICAgc3FydCA9IE1hdGguc3FydCxcbiAgICBjb3MgPSBNYXRoLmNvcyxcbiAgICBzaW4gPSBNYXRoLnNpblxuXG4gIENvbXBsZXhBcnJheS5wcm90b3R5cGUuRkZUID0gZnVuY3Rpb24oKSB7XG4gICAgcmV0dXJuIEZGVCh0aGlzLCBmYWxzZSk7XG4gIH1cblxuICBleHBvcnRzLkZGVCA9IGZ1bmN0aW9uKGlucHV0KSB7XG4gICAgcmV0dXJuIGVuc3VyZUNvbXBsZXhBcnJheShpbnB1dCkuRkZUKClcbiAgfVxuXG4gIENvbXBsZXhBcnJheS5wcm90b3R5cGUuSW52RkZUID0gZnVuY3Rpb24oKSB7XG4gICAgcmV0dXJuIEZGVCh0aGlzLCB0cnVlKVxuICB9XG5cbiAgZXhwb3J0cy5JbnZGRlQgPSBmdW5jdGlvbihpbnB1dCkge1xuICAgIHJldHVybiBlbnN1cmVDb21wbGV4QXJyYXkoaW5wdXQpLkludkZGVCgpXG4gIH1cblxuICAvLyBBcHBsaWVzIGEgZnJlcXVlbmN5LXNwYWNlIGZpbHRlciB0byBpbnB1dCwgYW5kIHJldHVybnMgdGhlIHJlYWwtc3BhY2VcbiAgLy8gZmlsdGVyZWQgaW5wdXQuXG4gIC8vIGZpbHRlcmVyIGFjY2VwdHMgZnJlcSwgaSwgbiBhbmQgbW9kaWZpZXMgZnJlcS5yZWFsIGFuZCBmcmVxLmltYWcuXG4gIENvbXBsZXhBcnJheS5wcm90b3R5cGUuZnJlcXVlbmN5TWFwID0gZnVuY3Rpb24oZmlsdGVyZXIpIHtcbiAgICByZXR1cm4gdGhpcy5GRlQoKS5tYXAoZmlsdGVyZXIpLkludkZGVCgpXG4gIH1cblxuICBleHBvcnRzLmZyZXF1ZW5jeU1hcCA9IGZ1bmN0aW9uKGlucHV0LCBmaWx0ZXJlcikge1xuICAgIHJldHVybiBlbnN1cmVDb21wbGV4QXJyYXkoaW5wdXQpLmZyZXF1ZW5jeU1hcChmaWx0ZXJlcilcbiAgfVxuXG4gIGZ1bmN0aW9uIGVuc3VyZUNvbXBsZXhBcnJheShpbnB1dCkge1xuICAgIHJldHVybiBjb21wbGV4X2FycmF5LmlzQ29tcGxleEFycmF5KGlucHV0KSAmJiBpbnB1dCB8fFxuICAgICAgICBuZXcgQ29tcGxleEFycmF5KGlucHV0KVxuICB9XG5cbiAgZnVuY3Rpb24gRkZUKGlucHV0LCBpbnZlcnNlKSB7XG4gICAgdmFyIG4gPSBpbnB1dC5sZW5ndGhcblxuICAgIGlmIChuICYgKG4gLSAxKSkge1xuICAgICAgcmV0dXJuIEZGVF9SZWN1cnNpdmUoaW5wdXQsIGludmVyc2UpXG4gICAgfSBlbHNlIHtcbiAgICAgIHJldHVybiBGRlRfMl9JdGVyYXRpdmUoaW5wdXQsIGludmVyc2UpXG4gICAgfVxuICB9XG5cbiAgZnVuY3Rpb24gRkZUX1JlY3Vyc2l2ZShpbnB1dCwgaW52ZXJzZSkge1xuICAgIHZhclxuICAgICAgbiA9IGlucHV0Lmxlbmd0aCxcbiAgICAgIC8vIENvdW50ZXJzLlxuICAgICAgaSwgaixcbiAgICAgIG91dHB1dCxcbiAgICAgIC8vIENvbXBsZXggbXVsdGlwbGllciBhbmQgaXRzIGRlbHRhLlxuICAgICAgZl9yLCBmX2ksIGRlbF9mX3IsIGRlbF9mX2ksXG4gICAgICAvLyBMb3dlc3QgZGl2aXNvciBhbmQgcmVtYWluZGVyLlxuICAgICAgcCwgbSxcbiAgICAgIG5vcm1hbGlzYXRpb24sXG4gICAgICByZWN1cnNpdmVfcmVzdWx0LFxuICAgICAgX3N3YXAsIF9yZWFsLCBfaW1hZ1xuXG4gICAgaWYgKG4gPT09IDEpIHtcbiAgICAgIHJldHVybiBpbnB1dFxuICAgIH1cblxuICAgIG91dHB1dCA9IG5ldyBDb21wbGV4QXJyYXkobiwgaW5wdXQuQXJyYXlUeXBlKVxuXG4gICAgLy8gVXNlIHRoZSBsb3dlc3Qgb2RkIGZhY3Rvciwgc28gd2UgYXJlIGFibGUgdG8gdXNlIEZGVF8yX0l0ZXJhdGl2ZSBpbiB0aGVcbiAgICAvLyByZWN1cnNpdmUgdHJhbnNmb3JtcyBvcHRpbWFsbHkuXG4gICAgcCA9IExvd2VzdE9kZEZhY3RvcihuKVxuICAgIG0gPSBuIC8gcFxuICAgIG5vcm1hbGlzYXRpb24gPSAxIC8gc3FydChwKVxuICAgIHJlY3Vyc2l2ZV9yZXN1bHQgPSBuZXcgQ29tcGxleEFycmF5KG0sIGlucHV0LkFycmF5VHlwZSlcblxuICAgIC8vIExvb3BzIGdvIGxpa2UgTyhuIM6jIHBfaSksIHdoZXJlIHBfaSBhcmUgdGhlIHByaW1lIGZhY3RvcnMgb2Ygbi5cbiAgICAvLyBmb3IgYSBwb3dlciBvZiBhIHByaW1lLCBwLCB0aGlzIHJlZHVjZXMgdG8gTyhuIHAgbG9nX3AgbilcbiAgICBmb3IoaiA9IDA7IGogPCBwOyBqKyspIHtcbiAgICAgIGZvcihpID0gMDsgaSA8IG07IGkrKykge1xuICAgICAgICByZWN1cnNpdmVfcmVzdWx0LnJlYWxbaV0gPSBpbnB1dC5yZWFsW2kgKiBwICsgal1cbiAgICAgICAgcmVjdXJzaXZlX3Jlc3VsdC5pbWFnW2ldID0gaW5wdXQuaW1hZ1tpICogcCArIGpdXG4gICAgICB9XG4gICAgICAvLyBEb24ndCBnbyBkZWVwZXIgdW5sZXNzIG5lY2Vzc2FyeSB0byBzYXZlIGFsbG9jcy5cbiAgICAgIGlmIChtID4gMSkge1xuICAgICAgICByZWN1cnNpdmVfcmVzdWx0ID0gRkZUKHJlY3Vyc2l2ZV9yZXN1bHQsIGludmVyc2UpXG4gICAgICB9XG5cbiAgICAgIGRlbF9mX3IgPSBjb3MoMipQSSpqL24pXG4gICAgICBkZWxfZl9pID0gKGludmVyc2UgPyAtMSA6IDEpICogc2luKDIqUEkqai9uKVxuICAgICAgZl9yID0gMVxuICAgICAgZl9pID0gMFxuXG4gICAgICBmb3IoaSA9IDA7IGkgPCBuOyBpKyspIHtcbiAgICAgICAgX3JlYWwgPSByZWN1cnNpdmVfcmVzdWx0LnJlYWxbaSAlIG1dXG4gICAgICAgIF9pbWFnID0gcmVjdXJzaXZlX3Jlc3VsdC5pbWFnW2kgJSBtXVxuXG4gICAgICAgIG91dHB1dC5yZWFsW2ldICs9IGZfciAqIF9yZWFsIC0gZl9pICogX2ltYWdcbiAgICAgICAgb3V0cHV0LmltYWdbaV0gKz0gZl9yICogX2ltYWcgKyBmX2kgKiBfcmVhbFxuXG4gICAgICAgIF9zd2FwID0gZl9yICogZGVsX2ZfciAtIGZfaSAqIGRlbF9mX2lcbiAgICAgICAgZl9pID0gZl9yICogZGVsX2ZfaSArIGZfaSAqIGRlbF9mX3JcbiAgICAgICAgZl9yID0gX3N3YXBcbiAgICAgIH1cbiAgICB9XG5cbiAgICAvLyBDb3B5IGJhY2sgdG8gaW5wdXQgdG8gbWF0Y2ggRkZUXzJfSXRlcmF0aXZlIGluLXBsYWNlbmVzc1xuICAgIC8vIFRPRE86IGZhc3RlciB3YXkgb2YgbWFraW5nIHRoaXMgaW4tcGxhY2U/XG4gICAgZm9yKGkgPSAwOyBpIDwgbjsgaSsrKSB7XG4gICAgICBpbnB1dC5yZWFsW2ldID0gbm9ybWFsaXNhdGlvbiAqIG91dHB1dC5yZWFsW2ldXG4gICAgICBpbnB1dC5pbWFnW2ldID0gbm9ybWFsaXNhdGlvbiAqIG91dHB1dC5pbWFnW2ldXG4gICAgfVxuXG4gICAgcmV0dXJuIGlucHV0XG4gIH1cblxuICBmdW5jdGlvbiBGRlRfMl9JdGVyYXRpdmUoaW5wdXQsIGludmVyc2UpIHtcbiAgICB2YXJcbiAgICAgIG4gPSBpbnB1dC5sZW5ndGgsXG4gICAgICAvLyBDb3VudGVycy5cbiAgICAgIGksIGosXG4gICAgICBvdXRwdXQsIG91dHB1dF9yLCBvdXRwdXRfaSxcbiAgICAgIC8vIENvbXBsZXggbXVsdGlwbGllciBhbmQgaXRzIGRlbHRhLlxuICAgICAgZl9yLCBmX2ksIGRlbF9mX3IsIGRlbF9mX2ksIHRlbXAsXG4gICAgICAvLyBUZW1wb3JhcnkgbG9vcCB2YXJpYWJsZXMuXG4gICAgICBsX2luZGV4LCByX2luZGV4LFxuICAgICAgbGVmdF9yLCBsZWZ0X2ksIHJpZ2h0X3IsIHJpZ2h0X2ksXG4gICAgICAvLyB3aWR0aCBvZiBlYWNoIHN1Yi1hcnJheSBmb3Igd2hpY2ggd2UncmUgaXRlcmF0aXZlbHkgY2FsY3VsYXRpbmcgRkZULlxuICAgICAgd2lkdGhcblxuICAgIG91dHB1dCA9IEJpdFJldmVyc2VDb21wbGV4QXJyYXkoaW5wdXQpXG4gICAgb3V0cHV0X3IgPSBvdXRwdXQucmVhbFxuICAgIG91dHB1dF9pID0gb3V0cHV0LmltYWdcbiAgICAvLyBMb29wcyBnbyBsaWtlIE8obiBsb2cgbik6XG4gICAgLy8gICB3aWR0aCB+IGxvZyBuOyBpLGogfiBuXG4gICAgd2lkdGggPSAxXG4gICAgd2hpbGUgKHdpZHRoIDwgbikge1xuICAgICAgZGVsX2ZfciA9IGNvcyhQSS93aWR0aClcbiAgICAgIGRlbF9mX2kgPSAoaW52ZXJzZSA/IC0xIDogMSkgKiBzaW4oUEkvd2lkdGgpXG4gICAgICBmb3IgKGkgPSAwOyBpIDwgbi8oMip3aWR0aCk7IGkrKykge1xuICAgICAgICBmX3IgPSAxXG4gICAgICAgIGZfaSA9IDBcbiAgICAgICAgZm9yIChqID0gMDsgaiA8IHdpZHRoOyBqKyspIHtcbiAgICAgICAgICBsX2luZGV4ID0gMippKndpZHRoICsgalxuICAgICAgICAgIHJfaW5kZXggPSBsX2luZGV4ICsgd2lkdGhcblxuICAgICAgICAgIGxlZnRfciA9IG91dHB1dF9yW2xfaW5kZXhdXG4gICAgICAgICAgbGVmdF9pID0gb3V0cHV0X2lbbF9pbmRleF1cbiAgICAgICAgICByaWdodF9yID0gZl9yICogb3V0cHV0X3Jbcl9pbmRleF0gLSBmX2kgKiBvdXRwdXRfaVtyX2luZGV4XVxuICAgICAgICAgIHJpZ2h0X2kgPSBmX2kgKiBvdXRwdXRfcltyX2luZGV4XSArIGZfciAqIG91dHB1dF9pW3JfaW5kZXhdXG5cbiAgICAgICAgICBvdXRwdXRfcltsX2luZGV4XSA9IFNRUlQxXzIgKiAobGVmdF9yICsgcmlnaHRfcilcbiAgICAgICAgICBvdXRwdXRfaVtsX2luZGV4XSA9IFNRUlQxXzIgKiAobGVmdF9pICsgcmlnaHRfaSlcbiAgICAgICAgICBvdXRwdXRfcltyX2luZGV4XSA9IFNRUlQxXzIgKiAobGVmdF9yIC0gcmlnaHRfcilcbiAgICAgICAgICBvdXRwdXRfaVtyX2luZGV4XSA9IFNRUlQxXzIgKiAobGVmdF9pIC0gcmlnaHRfaSlcbiAgICAgICAgICB0ZW1wID0gZl9yICogZGVsX2ZfciAtIGZfaSAqIGRlbF9mX2lcbiAgICAgICAgICBmX2kgPSBmX3IgKiBkZWxfZl9pICsgZl9pICogZGVsX2ZfclxuICAgICAgICAgIGZfciA9IHRlbXBcbiAgICAgICAgfVxuICAgICAgfVxuICAgICAgd2lkdGggPDw9IDFcbiAgICB9XG5cbiAgICByZXR1cm4gb3V0cHV0XG4gIH1cblxuICBmdW5jdGlvbiBCaXRSZXZlcnNlSW5kZXgoaW5kZXgsIG4pIHtcbiAgICB2YXIgYml0cmV2ZXJzZWRfaW5kZXggPSAwXG5cbiAgICB3aGlsZSAobiA+IDEpIHtcbiAgICAgIGJpdHJldmVyc2VkX2luZGV4IDw8PSAxXG4gICAgICBiaXRyZXZlcnNlZF9pbmRleCArPSBpbmRleCAmIDFcbiAgICAgIGluZGV4ID4+PSAxXG4gICAgICBuID4+PSAxXG4gICAgfVxuICAgIHJldHVybiBiaXRyZXZlcnNlZF9pbmRleFxuICB9XG5cbiAgZnVuY3Rpb24gQml0UmV2ZXJzZUNvbXBsZXhBcnJheShhcnJheSkge1xuICAgIHZhciBuID0gYXJyYXkubGVuZ3RoLFxuICAgICAgICBmbGlwcyA9IHt9LFxuICAgICAgICBzd2FwLFxuICAgICAgICBpXG5cbiAgICBmb3IoaSA9IDA7IGkgPCBuOyBpKyspIHtcbiAgICAgIHZhciByX2kgPSBCaXRSZXZlcnNlSW5kZXgoaSwgbilcblxuICAgICAgaWYgKGZsaXBzLmhhc093blByb3BlcnR5KGkpIHx8IGZsaXBzLmhhc093blByb3BlcnR5KHJfaSkpIGNvbnRpbnVlXG5cbiAgICAgIHN3YXAgPSBhcnJheS5yZWFsW3JfaV1cbiAgICAgIGFycmF5LnJlYWxbcl9pXSA9IGFycmF5LnJlYWxbaV1cbiAgICAgIGFycmF5LnJlYWxbaV0gPSBzd2FwXG5cbiAgICAgIHN3YXAgPSBhcnJheS5pbWFnW3JfaV1cbiAgICAgIGFycmF5LmltYWdbcl9pXSA9IGFycmF5LmltYWdbaV1cbiAgICAgIGFycmF5LmltYWdbaV0gPSBzd2FwXG5cbiAgICAgIGZsaXBzW2ldID0gZmxpcHNbcl9pXSA9IHRydWVcbiAgICB9XG5cbiAgICByZXR1cm4gYXJyYXlcbiAgfVxuXG4gIGZ1bmN0aW9uIExvd2VzdE9kZEZhY3RvcihuKSB7XG4gICAgdmFyIGZhY3RvciA9IDMsXG4gICAgICAgIHNxcnRfbiA9IHNxcnQobilcblxuICAgIHdoaWxlKGZhY3RvciA8PSBzcXJ0X24pIHtcbiAgICAgIGlmIChuICUgZmFjdG9yID09PSAwKSByZXR1cm4gZmFjdG9yXG4gICAgICBmYWN0b3IgPSBmYWN0b3IgKyAyXG4gICAgfVxuICAgIHJldHVybiBuXG4gIH1cblxufShcbiAgdHlwZW9mIGV4cG9ydHMgPT09ICd1bmRlZmluZWQnICYmICh0aGlzLmZmdCA9IHt9KSB8fCBleHBvcnRzLFxuICB0eXBlb2YgcmVxdWlyZSA9PT0gJ3VuZGVmaW5lZCcgJiYgKHRoaXMuY29tcGxleF9hcnJheSkgfHxcbiAgICByZXF1aXJlKCcuL2NvbXBsZXhfYXJyYXknKVxuKVxuIiwibW9kdWxlLmV4cG9ydHMgPSBmdW5jdGlvbihidWZmZXJTaXplLCBtKXtcbiAgcmV0dXJuIG0uYW1wU3BlY3RydW07XG59OyIsIm1vZHVsZS5leHBvcnRzID0gZnVuY3Rpb24oYnVmZmVyU2l6ZSwgbSkge1xuICByZXR1cm4gbS5zaWduYWw7XG59OyIsIm1vZHVsZS5leHBvcnRzID0gZnVuY3Rpb24oYnVmZmVyU2l6ZSwgbSkge1xuICByZXR1cm4gbS5jb21wbGV4U3BlY3RydW07XG59OyIsIm1vZHVsZS5leHBvcnRzID0gZnVuY3Rpb24oYnVmZmVyU2l6ZSwgbSkge1xuICB2YXIgZW5lcmd5ID0gMDtcbiAgZm9yICh2YXIgaSA9IDA7IGkgPCBtLnNpZ25hbC5sZW5ndGg7IGkrKykge1xuICAgIGVuZXJneSArPSBNYXRoLnBvdyhNYXRoLmFicyhtLnNpZ25hbFtpXSksIDIpO1xuICB9XG4gIHJldHVybiBlbmVyZ3k7XG59OyIsIm1vZHVsZS5leHBvcnRzID0ge1xuICBidWZmZXI6IHJlcXVpcmUoJy4vYnVmZmVyJyksXG4gIHJtczogcmVxdWlyZSgnLi9ybXMnKSxcbiAgZW5lcmd5OiByZXF1aXJlKCcuL2VuZXJneScpLFxuICBjb21wbGV4U3BlY3RydW06IHJlcXVpcmUoJy4vY29tcGxleFNwZWN0cnVtJyksXG4gIHNwZWN0cmFsU2xvcGU6IHJlcXVpcmUoJy4vc3BlY3RyYWxTbG9wZScpLFxuICBzcGVjdHJhbENlbnRyb2lkOiByZXF1aXJlKCcuL3NwZWN0cmFsQ2VudHJvaWQnKSxcbiAgc3BlY3RyYWxSb2xsb2ZmOiByZXF1aXJlKCcuL3NwZWN0cmFsUm9sbG9mZicpLFxuICBzcGVjdHJhbEZsYXRuZXNzOiByZXF1aXJlKCcuL3NwZWN0cmFsRmxhdG5lc3MnKSxcbiAgc3BlY3RyYWxTcHJlYWQ6IHJlcXVpcmUoJy4vc3BlY3RyYWxTcHJlYWQnKSxcbiAgc3BlY3RyYWxTa2V3bmVzczogcmVxdWlyZSgnLi9zcGVjdHJhbFNrZXduZXNzJyksXG4gIHNwZWN0cmFsS3VydG9zaXM6IHJlcXVpcmUoJy4vc3BlY3RyYWxLdXJ0b3NpcycpLFxuICBhbXBsaXR1ZGVTcGVjdHJ1bTogcmVxdWlyZSgnLi9hbXBsaXR1ZGVTcGVjdHJ1bScpLFxuICB6Y3I6IHJlcXVpcmUoJy4vemNyJyksXG4gIHBvd2VyU3BlY3RydW06IHJlcXVpcmUoJy4vcG93ZXJTcGVjdHJ1bScpLFxuICBsb3VkbmVzczogcmVxdWlyZSgnLi9sb3VkbmVzcycpLFxuICBwZXJjZXB0dWFsU3ByZWFkOiByZXF1aXJlKCcuL3BlcmNlcHR1YWxTcHJlYWQnKSxcbiAgcGVyY2VwdHVhbFNoYXJwbmVzczogcmVxdWlyZSgnLi9wZXJjZXB0dWFsU2hhcnBuZXNzJyksXG4gIG1mY2M6IHJlcXVpcmUoJy4vbWZjYycpXG59OyIsIm1vZHVsZS5leHBvcnRzID0gZnVuY3Rpb24oYnVmZmVyU2l6ZSwgbSkge1xuICB2YXIgYmFya1NjYWxlID0gbmV3IEZsb2F0MzJBcnJheShtLmFtcFNwZWN0cnVtLmxlbmd0aCk7XG4gIHZhciBOVU1fQkFSS19CQU5EUyA9IDI0O1xuICB2YXIgc3BlY2lmaWMgPSBuZXcgRmxvYXQzMkFycmF5KE5VTV9CQVJLX0JBTkRTKTtcbiAgdmFyIHRvdCA9IDA7XG4gIHZhciBub3JtYWxpc2VkU3BlY3RydW0gPSBtLmFtcFNwZWN0cnVtO1xuICB2YXIgYmJMaW1pdHMgPSBuZXcgSW50MzJBcnJheShOVU1fQkFSS19CQU5EUyArIDEpO1xuXG4gIGZvciAodmFyIGkgPSAwOyBpIDwgYmFya1NjYWxlLmxlbmd0aDsgaSsrKSB7XG4gICAgYmFya1NjYWxlW2ldID0gaSAqIG0uYXVkaW9Db250ZXh0LnNhbXBsZVJhdGUgLyAoYnVmZmVyU2l6ZSk7XG4gICAgYmFya1NjYWxlW2ldID0gMTMgKiBNYXRoLmF0YW4oYmFya1NjYWxlW2ldIC8gMTMxNS44KSArIDMuNSAqIE1hdGguYXRhbihNYXRoLnBvdygoYmFya1NjYWxlW2ldIC8gNzUxOCksIDIpKTtcbiAgfVxuXG5cbiAgYmJMaW1pdHNbMF0gPSAwO1xuICB2YXIgY3VycmVudEJhbmRFbmQgPSBiYXJrU2NhbGVbbS5hbXBTcGVjdHJ1bS5sZW5ndGggLSAxXSAvIE5VTV9CQVJLX0JBTkRTO1xuICB2YXIgY3VycmVudEJhbmQgPSAxO1xuICBmb3IgKHZhciBpID0gMDsgaSA8IG0uYW1wU3BlY3RydW0ubGVuZ3RoOyBpKyspIHtcbiAgICB3aGlsZSAoYmFya1NjYWxlW2ldID4gY3VycmVudEJhbmRFbmQpIHtcbiAgICAgIGJiTGltaXRzW2N1cnJlbnRCYW5kKytdID0gaTtcbiAgICAgIGN1cnJlbnRCYW5kRW5kID0gY3VycmVudEJhbmQgKiBiYXJrU2NhbGVbbS5hbXBTcGVjdHJ1bS5sZW5ndGggLSAxXSAvIE5VTV9CQVJLX0JBTkRTO1xuICAgIH1cbiAgfVxuXG4gIGJiTGltaXRzW05VTV9CQVJLX0JBTkRTXSA9IG0uYW1wU3BlY3RydW0ubGVuZ3RoIC0gMTtcblxuICAvL3Byb2Nlc3NcblxuICBmb3IgKHZhciBpID0gMDsgaSA8IE5VTV9CQVJLX0JBTkRTOyBpKyspIHtcbiAgICB2YXIgc3VtID0gMDtcbiAgICBmb3IgKHZhciBqID0gYmJMaW1pdHNbaV07IGogPCBiYkxpbWl0c1tpICsgMV07IGorKykge1xuXG4gICAgICBzdW0gKz0gbm9ybWFsaXNlZFNwZWN0cnVtW2pdO1xuICAgIH1cbiAgICBzcGVjaWZpY1tpXSA9IE1hdGgucG93KHN1bSwgMC4yMyk7XG4gIH1cblxuICAvL2dldCB0b3RhbCBsb3VkbmVzc1xuICBmb3IgKHZhciBpID0gMDsgaSA8IHNwZWNpZmljLmxlbmd0aDsgaSsrKSB7XG4gICAgdG90ICs9IHNwZWNpZmljW2ldO1xuICB9XG4gIHJldHVybiB7XG4gICAgXCJzcGVjaWZpY1wiOiBzcGVjaWZpYyxcbiAgICBcInRvdGFsXCI6IHRvdFxuICB9O1xufTsiLCIvL3VzZWQgdHV0b3JpYWwgZnJvbSBodHRwOi8vcHJhY3RpY2FsY3J5cHRvZ3JhcGh5LmNvbS9taXNjZWxsYW5lb3VzL21hY2hpbmUtbGVhcm5pbmcvZ3VpZGUtbWVsLWZyZXF1ZW5jeS1jZXBzdHJhbC1jb2VmZmljaWVudHMtbWZjY3MvXG5cbnZhciBwb3dlclNwZWN0cnVtID0gcmVxdWlyZSgnLi9wb3dlclNwZWN0cnVtJyk7XG5cbm1vZHVsZS5leHBvcnRzID0gZnVuY3Rpb24oYnVmZmVyU2l6ZSwgbSkge1xuICB2YXIgcG93U3BlYyA9IHBvd2VyU3BlY3RydW0oYnVmZmVyU2l6ZSwgbSk7XG4gIHZhciBmcmVxVG9NZWwgPSBmdW5jdGlvbihmcmVxVmFsdWUpIHtcbiAgICB2YXIgbWVsVmFsdWUgPSAxMTI1ICogTWF0aC5sb2coMSArIChmcmVxVmFsdWUgLyA3MDApKTtcbiAgICByZXR1cm4gbWVsVmFsdWU7XG4gIH07XG4gIHZhciBtZWxUb0ZyZXEgPSBmdW5jdGlvbihtZWxWYWx1ZSkge1xuICAgIHZhciBmcmVxVmFsdWUgPSA3MDAgKiAoTWF0aC5leHAobWVsVmFsdWUgLyAxMTI1KSAtIDEpO1xuICAgIHJldHVybiBmcmVxVmFsdWU7XG4gIH07XG4gIHZhciBudW1GaWx0ZXJzID0gMjY7IC8vMjYgZmlsdGVycyBpcyBzdGFuZGFyZFxuICB2YXIgbWVsVmFsdWVzID0gbmV3IEZsb2F0MzJBcnJheShudW1GaWx0ZXJzICsgMik7IC8vdGhlICsyIGlzIHRoZSB1cHBlciBhbmQgbG93ZXIgbGltaXRzXG4gIHZhciBtZWxWYWx1ZXNJbkZyZXEgPSBuZXcgRmxvYXQzMkFycmF5KG51bUZpbHRlcnMgKyAyKTtcbiAgLy9HZW5lcmF0ZSBsaW1pdHMgaW4gSHogLSBmcm9tIDAgdG8gdGhlIG55cXVpc3QuXG4gIHZhciBsb3dlckxpbWl0RnJlcSA9IDA7XG4gIHZhciB1cHBlckxpbWl0RnJlcSA9IGF1ZGlvQ29udGV4dC5zYW1wbGVSYXRlIC8gMjtcbiAgLy9Db252ZXJ0IHRoZSBsaW1pdHMgdG8gTWVsXG4gIHZhciBsb3dlckxpbWl0TWVsID0gZnJlcVRvTWVsKGxvd2VyTGltaXRGcmVxKTtcbiAgdmFyIHVwcGVyTGltaXRNZWwgPSBmcmVxVG9NZWwodXBwZXJMaW1pdEZyZXEpO1xuICAvL0ZpbmQgdGhlIHJhbmdlXG4gIHZhciByYW5nZSA9IHVwcGVyTGltaXRNZWwgLSBsb3dlckxpbWl0TWVsO1xuICAvL0ZpbmQgdGhlIHJhbmdlIGFzIHBhcnQgb2YgdGhlIGxpbmVhciBpbnRlcnBvbGF0aW9uXG4gIHZhciB2YWx1ZVRvQWRkID0gcmFuZ2UgLyAobnVtRmlsdGVycyArIDEpO1xuXG4gIHZhciBmZnRCaW5zT2ZGcmVxID0gQXJyYXkobnVtRmlsdGVycyArIDIpO1xuXG4gIGZvciAodmFyIGkgPSAwOyBpIDwgbWVsVmFsdWVzLmxlbmd0aDsgaSsrKSB7XG4gICAgLy9Jbml0aWFsaXNpbmcgdGhlIG1lbCBmcmVxdWVuY2llcyAtIHRoZXkgYXJlIGp1c3QgYSBsaW5lYXIgaW50ZXJwb2xhdGlvbiBiZXR3ZWVuIHRoZSBsb3dlciBhbmQgdXBwZXIgbGltaXRzLlxuICAgIG1lbFZhbHVlc1tpXSA9IGkgKiB2YWx1ZVRvQWRkO1xuICAgIC8vQ29udmVydCBiYWNrIHRvIEh6XG4gICAgbWVsVmFsdWVzSW5GcmVxW2ldID0gbWVsVG9GcmVxKG1lbFZhbHVlc1tpXSk7XG4gICAgLy9GaW5kIHRoZSBjb3JyZXNwb25kaW5nIGJpbnNcbiAgICBmZnRCaW5zT2ZGcmVxW2ldID0gTWF0aC5mbG9vcigoYnVmZmVyU2l6ZSArIDEpICogbWVsVmFsdWVzSW5GcmVxW2ldIC8gYXVkaW9Db250ZXh0LnNhbXBsZVJhdGUpO1xuICB9XG5cbiAgdmFyIGZpbHRlckJhbmsgPSBBcnJheShudW1GaWx0ZXJzKTtcbiAgZm9yICh2YXIgaiA9IDA7IGogPCBmaWx0ZXJCYW5rLmxlbmd0aDsgaisrKSB7XG4gICAgLy9jcmVhdGluZyBhIHR3byBkaW1lbnNpb25hbCBhcnJheSBvZiBzaXplIG51bUZpbHRlcyAqIChidWZmZXJzaXplLzIpKzEgYW5kIHByZS1wb3B1bGF0aW5nIHRoZSBhcnJheXMgd2l0aCAwcy5cbiAgICBmaWx0ZXJCYW5rW2pdID0gQXJyYXkuYXBwbHkobnVsbCwgbmV3IEFycmF5KChidWZmZXJTaXplIC8gMikgKyAxKSkubWFwKE51bWJlci5wcm90b3R5cGUudmFsdWVPZiwgMCk7XG4gICAgLy9jcmVhdGluZyB0aGUgbG93ZXIgYW5kIHVwcGVyIHNsb3BlcyBmb3IgZWFjaCBiaW5cbiAgICBmb3IgKHZhciBpID0gZmZ0Qmluc09mRnJlcVtqXTsgaSA8IGZmdEJpbnNPZkZyZXFbaiArIDFdOyBpKyspIHtcbiAgICAgIGZpbHRlckJhbmtbal1baV0gPSAoaSAtIGZmdEJpbnNPZkZyZXFbal0pIC8gKGZmdEJpbnNPZkZyZXFbaiArIDFdIC0gZmZ0Qmluc09mRnJlcVtqXSk7XG4gICAgfVxuICAgIGZvciAodmFyIGkgPSBmZnRCaW5zT2ZGcmVxW2ogKyAxXTsgaSA8IGZmdEJpbnNPZkZyZXFbaiArIDJdOyBpKyspIHtcbiAgICAgIGZpbHRlckJhbmtbal1baV0gPSAoZmZ0Qmluc09mRnJlcVtqICsgMl0gLSBpKSAvIChmZnRCaW5zT2ZGcmVxW2ogKyAyXSAtIGZmdEJpbnNPZkZyZXFbaiArIDFdKTtcbiAgICB9XG4gIH1cblxuICB2YXIgbG9nZ2VkTWVsQmFuZHMgPSBuZXcgRmxvYXQzMkFycmF5KG51bUZpbHRlcnMpO1xuICBmb3IgKHZhciBpID0gMDsgaSA8IGxvZ2dlZE1lbEJhbmRzLmxlbmd0aDsgaSsrKSB7XG4gICAgbG9nZ2VkTWVsQmFuZHNbaV0gPSAwO1xuICAgIGZvciAodmFyIGogPSAwOyBqIDwgKGJ1ZmZlclNpemUgLyAyKTsgaisrKSB7XG4gICAgICAvL3BvaW50IG11bHRpcGxpY2F0aW9uIGJldHdlZW4gcG93ZXIgc3BlY3RydW0gYW5kIGZpbHRlcmJhbmtzLlxuICAgICAgZmlsdGVyQmFua1tpXVtqXSA9IGZpbHRlckJhbmtbaV1bal0gKiBwb3dTcGVjW2pdO1xuXG4gICAgICAvL3N1bW1pbmcgdXAgYWxsIG9mIHRoZSBjb2VmZmljaWVudHMgaW50byBvbmUgYXJyYXlcbiAgICAgIGxvZ2dlZE1lbEJhbmRzW2ldICs9IGZpbHRlckJhbmtbaV1bal07XG4gICAgfVxuICAgIC8vbG9nIGVhY2ggY29lZmZpY2llbnRcbiAgICBsb2dnZWRNZWxCYW5kc1tpXSA9IE1hdGgubG9nKGxvZ2dlZE1lbEJhbmRzW2ldKTtcbiAgfVxuXG4gIC8vZGN0XG4gIHZhciBrID0gTWF0aC5QSSAvIG51bUZpbHRlcnM7XG4gIHZhciB3MSA9IDEuMCAvIE1hdGguc3FydChudW1GaWx0ZXJzKTtcbiAgdmFyIHcyID0gTWF0aC5zcXJ0KDIuMCAvIG51bUZpbHRlcnMpO1xuICB2YXIgbnVtQ29lZmZzID0gMTM7XG4gIHZhciBkY3RNYXRyaXggPSBuZXcgRmxvYXQzMkFycmF5KG51bUNvZWZmcyAqIG51bUZpbHRlcnMpO1xuXG4gIGZvciAodmFyIGkgPSAwOyBpIDwgbnVtQ29lZmZzOyBpKyspIHtcbiAgICBmb3IgKHZhciBqID0gMDsgaiA8IG51bUZpbHRlcnM7IGorKykge1xuICAgICAgdmFyIGlkeCA9IGkgKyAoaiAqIG51bUNvZWZmcyk7XG4gICAgICBpZiAoaSA9PSAwKSB7XG4gICAgICAgIGRjdE1hdHJpeFtpZHhdID0gdzEgKiBNYXRoLmNvcyhrICogKGkgKyAxKSAqIChqICsgMC41KSk7XG4gICAgICB9IGVsc2Uge1xuICAgICAgICBkY3RNYXRyaXhbaWR4XSA9IHcyICogTWF0aC5jb3MoayAqIChpICsgMSkgKiAoaiArIDAuNSkpO1xuICAgICAgfVxuICAgIH1cbiAgfVxuXG4gIHZhciBtZmNjcyA9IG5ldyBGbG9hdDMyQXJyYXkobnVtQ29lZmZzKTtcbiAgZm9yICh2YXIgayA9IDA7IGsgPCBudW1Db2VmZnM7IGsrKykge1xuICAgIHZhciB2ID0gMDtcbiAgICBmb3IgKHZhciBuID0gMDsgbiA8IG51bUZpbHRlcnM7IG4rKykge1xuICAgICAgdmFyIGlkeCA9IGsgKyAobiAqIG51bUNvZWZmcyk7XG4gICAgICB2ICs9IChkY3RNYXRyaXhbaWR4XSAqIGxvZ2dlZE1lbEJhbmRzW25dKTtcbiAgICB9XG4gICAgbWZjY3Nba10gPSB2IC8gbnVtQ29lZmZzO1xuICB9XG5cbiAgcmV0dXJuIG1mY2NzO1xufTsiLCJtb2R1bGUuZXhwb3J0cyA9IGZ1bmN0aW9uKGJ1ZmZlclNpemUsIG0pIHtcbiAgdmFyIGxvdWRuZXNzID0gbS5mZWF0dXJlRXh0cmFjdG9ycy5sb3VkbmVzcyhidWZmZXJTaXplLCBtKTtcbiAgdmFyIHNwZWMgPSBsb3VkbmVzcy5zcGVjaWZpYztcbiAgdmFyIG91dHB1dCA9IDA7XG5cbiAgZm9yICh2YXIgaSA9IDA7IGkgPCBzcGVjLmxlbmd0aDsgaSsrKSB7XG4gICAgaWYgKGkgPCAxNSkge1xuICAgICAgb3V0cHV0ICs9IChpICsgMSkgKiBzcGVjW2kgKyAxXTtcbiAgICB9IGVsc2Uge1xuICAgICAgb3V0cHV0ICs9IDAuMDY2ICogTWF0aC5leHAoMC4xNzEgKiAoaSArIDEpKTtcbiAgICB9XG4gIH1cbiAgb3V0cHV0ICo9IDAuMTEgLyBsb3VkbmVzcy50b3RhbDtcblxuICByZXR1cm4gb3V0cHV0O1xufTsiLCJtb2R1bGUuZXhwb3J0cyA9IGZ1bmN0aW9uKGJ1ZmZlclNpemUsIG0pIHtcbiAgdmFyIGxvdWRuZXNzID0gbS5mZWF0dXJlRXh0cmFjdG9ycy5sb3VkbmVzcyhidWZmZXJTaXplLCBtKTtcblxuICB2YXIgbWF4ID0gMDtcbiAgZm9yICh2YXIgaSA9IDA7IGkgPCBsb3VkbmVzcy5zcGVjaWZpYy5sZW5ndGg7IGkrKykge1xuICAgIGlmIChsb3VkbmVzcy5zcGVjaWZpY1tpXSA+IG1heCkge1xuICAgICAgbWF4ID0gbG91ZG5lc3Muc3BlY2lmaWNbaV07XG4gICAgfVxuICB9XG5cbiAgdmFyIHNwcmVhZCA9IE1hdGgucG93KChsb3VkbmVzcy50b3RhbCAtIG1heCkgLyBsb3VkbmVzcy50b3RhbCwgMik7XG5cbiAgcmV0dXJuIHNwcmVhZDtcbn07IiwibW9kdWxlLmV4cG9ydHMgPSBmdW5jdGlvbihidWZmZXJTaXplLCBtKSB7XG4gIHZhciBwb3dlclNwZWN0cnVtID0gbmV3IEZsb2F0MzJBcnJheShtLmFtcFNwZWN0cnVtLmxlbmd0aCk7XG4gIGZvciAodmFyIGkgPSAwOyBpIDwgcG93ZXJTcGVjdHJ1bS5sZW5ndGg7IGkrKykge1xuICAgIHBvd2VyU3BlY3RydW1baV0gPSBNYXRoLnBvdyhtLmFtcFNwZWN0cnVtW2ldLCAyKTtcbiAgfVxuICByZXR1cm4gcG93ZXJTcGVjdHJ1bTtcbn07IiwibW9kdWxlLmV4cG9ydHMgPSBmdW5jdGlvbihidWZmZXJTaXplLCBtKSB7XG5cbiAgdmFyIHJtcyA9IDA7XG4gIGZvciAodmFyIGkgPSAwOyBpIDwgbS5zaWduYWwubGVuZ3RoOyBpKyspIHtcbiAgICBybXMgKz0gTWF0aC5wb3cobS5zaWduYWxbaV0sIDIpO1xuICB9XG4gIHJtcyA9IHJtcyAvIG0uc2lnbmFsLmxlbmd0aDtcbiAgcm1zID0gTWF0aC5zcXJ0KHJtcyk7XG5cbiAgcmV0dXJuIHJtcztcbn07XG4iLCJtb2R1bGUuZXhwb3J0cyA9IGZ1bmN0aW9uKGJ1ZmZlclNpemUsIG0pIHtcbiAgcmV0dXJuIMK1KDEsIG0uYW1wU3BlY3RydW0pO1xufTsiLCJtb2R1bGUuZXhwb3J0cyA9IGZ1bmN0aW9uKGJ1ZmZlclNpemUsIG0pIHtcbiAgdmFyIGFtcHNwZWMgPSBtLmFtcFNwZWN0cnVtO1xuICB2YXIgbnVtZXJhdG9yID0gMDtcbiAgdmFyIGRlbm9taW5hdG9yID0gMDtcbiAgZm9yICh2YXIgaSA9IDA7IGkgPCBhbXBzcGVjLmxlbmd0aDsgaSsrKSB7XG4gICAgbnVtZXJhdG9yICs9IE1hdGgubG9nKGFtcHNwZWNbaV0pO1xuICAgIGRlbm9taW5hdG9yICs9IGFtcHNwZWNbaV07XG4gIH1cbiAgcmV0dXJuIE1hdGguZXhwKG51bWVyYXRvciAvIGFtcHNwZWMubGVuZ3RoKSAqIGFtcHNwZWMubGVuZ3RoIC8gZGVub21pbmF0b3I7XG59OyIsIm1vZHVsZS5leHBvcnRzID0gZnVuY3Rpb24oYnVmZmVyU2l6ZSwgbSkge1xuICB2YXIgYW1wc3BlYyA9IG0uYW1wU3BlY3RydW07XG4gIHZhciDCtTEgPSDCtSgxLCBhbXBzcGVjKTtcbiAgdmFyIMK1MiA9IMK1KDIsIGFtcHNwZWMpO1xuICB2YXIgwrUzID0gwrUoMywgYW1wc3BlYyk7XG4gIHZhciDCtTQgPSDCtSg0LCBhbXBzcGVjKTtcbiAgdmFyIG51bWVyYXRvciA9IC0zICogTWF0aC5wb3cowrUxLCA0KSArIDYgKiDCtTEgKiDCtTIgLSA0ICogwrUxICogwrUzICsgwrU0O1xuICB2YXIgZGVub21pbmF0b3IgPSBNYXRoLnBvdyhNYXRoLnNxcnQowrUyIC0gTWF0aC5wb3cowrUxLCAyKSksIDQpO1xuICByZXR1cm4gbnVtZXJhdG9yIC8gZGVub21pbmF0b3I7XG59OyIsIm1vZHVsZS5leHBvcnRzID0gZnVuY3Rpb24oYnVmZmVyU2l6ZSwgbSkge1xuICB2YXIgYW1wc3BlYyA9IG0uYW1wU3BlY3RydW07XG4gIC8vY2FsY3VsYXRlIG55cXVpc3QgYmluXG4gIHZhciBueXFCaW4gPSBtLmF1ZGlvQ29udGV4dC5zYW1wbGVSYXRlIC8gKDIgKiAoYW1wc3BlYy5sZW5ndGggLSAxKSk7XG4gIHZhciBlYyA9IDA7XG4gIGZvciAodmFyIGkgPSAwOyBpIDwgYW1wc3BlYy5sZW5ndGg7IGkrKykge1xuICAgIGVjICs9IGFtcHNwZWNbaV07XG4gIH1cbiAgdmFyIHRocmVzaG9sZCA9IDAuOTkgKiBlYztcbiAgdmFyIG4gPSBhbXBzcGVjLmxlbmd0aCAtIDE7XG4gIHdoaWxlIChlYyA+IHRocmVzaG9sZCAmJiBuID49IDApIHtcbiAgICBlYyAtPSBhbXBzcGVjW25dO1xuICAgIC0tbjtcbiAgfVxuICByZXR1cm4gKG4gKyAxKSAqIG55cUJpbjtcbn07IiwibW9kdWxlLmV4cG9ydHMgPSBmdW5jdGlvbihidWZmZXJTaXplLCBtLCBzcGVjdHJ1bSkge1xuICB2YXIgYW1wc3BlYyA9IG0uYW1wU3BlY3RydW07XG4gIHZhciDCtTEgPSDCtSgxLCBhbXBzcGVjKTtcbiAgdmFyIMK1MiA9IMK1KDIsIGFtcHNwZWMpO1xuICB2YXIgwrUzID0gwrUoMywgYW1wc3BlYyk7XG4gIHZhciBudW1lcmF0b3IgPSAyICogTWF0aC5wb3cowrUxLCAzKSAtIDMgKiDCtTEgKiDCtTIgKyDCtTM7XG4gIHZhciBkZW5vbWluYXRvciA9IE1hdGgucG93KE1hdGguc3FydCjCtTIgLSBNYXRoLnBvdyjCtTEsIDIpKSwgMyk7XG4gIHJldHVybiBudW1lcmF0b3IgLyBkZW5vbWluYXRvcjtcbn07IiwibW9kdWxlLmV4cG9ydHMgPSBmdW5jdGlvbihidWZmZXJTaXplLCBtKSB7XG4gIC8vbGluZWFyIHJlZ3Jlc3Npb25cbiAgdmFyIGFtcFN1bSA9IDA7XG4gIHZhciBmcmVxU3VtID0gMDtcbiAgdmFyIGZyZXFzID0gbmV3IEZsb2F0MzJBcnJheShtLmFtcFNwZWN0cnVtLmxlbmd0aCk7XG4gIHZhciBwb3dGcmVxU3VtID0gMDtcbiAgdmFyIGFtcEZyZXFTdW0gPSAwO1xuXG4gIGZvciAodmFyIGkgPSAwOyBpIDwgbS5hbXBTcGVjdHJ1bS5sZW5ndGg7IGkrKykge1xuICAgIGFtcFN1bSArPSBtLmFtcFNwZWN0cnVtW2ldO1xuICAgIHZhciBjdXJGcmVxID0gaSAqIG0uYXVkaW9Db250ZXh0LnNhbXBsZVJhdGUgLyBidWZmZXJTaXplO1xuICAgIGZyZXFzW2ldID0gY3VyRnJlcTtcbiAgICBwb3dGcmVxU3VtICs9IGN1ckZyZXEgKiBjdXJGcmVxO1xuICAgIGZyZXFTdW0gKz0gY3VyRnJlcTtcbiAgICBhbXBGcmVxU3VtICs9IGN1ckZyZXEgKiBtLmFtcFNwZWN0cnVtW2ldO1xuICB9XG4gIHJldHVybiAobS5hbXBTcGVjdHJ1bS5sZW5ndGggKiBhbXBGcmVxU3VtIC0gZnJlcVN1bSAqIGFtcFN1bSkgLyAoYW1wU3VtICogKHBvd0ZyZXFTdW0gLSBNYXRoLnBvdyhmcmVxU3VtLCAyKSkpO1xufTsiLCJtb2R1bGUuZXhwb3J0cyA9IGZ1bmN0aW9uKGJ1ZmZlclNpemUsIG0pIHtcbiAgdmFyIGFtcHNwZWMgPSBtLmFtcFNwZWN0cnVtO1xuICByZXR1cm4gTWF0aC5zcXJ0KMK1KDIsIGFtcHNwZWMpIC0gTWF0aC5wb3cowrUoMSwgYW1wc3BlYyksIDIpKTtcbn07IiwibW9kdWxlLmV4cG9ydHMgPSBmdW5jdGlvbihidWZmZXJTaXplLCBtKSB7XG4gIHZhciB6Y3IgPSAwO1xuICBmb3IgKHZhciBpID0gMDsgaSA8IG0uc2lnbmFsLmxlbmd0aDsgaSsrKSB7XG4gICAgaWYgKChtLnNpZ25hbFtpXSA+PSAwICYmIG0uc2lnbmFsW2kgKyAxXSA8IDApIHx8IChtLnNpZ25hbFtpXSA8IDAgJiYgbS5zaWduYWxbaSArIDFdID49IDApKSB7XG4gICAgICB6Y3IrKztcbiAgICB9XG4gIH1cbiAgcmV0dXJuIHpjcjtcbn07IiwiLy8gTWV5ZGEgSmF2YXNjcmlwdCBEU1AgbGlicmFyeVxuXG52YXIgQ29tcGxleEFycmF5ID0gcmVxdWlyZSgnLi4vbGliL2pzZmZ0L2NvbXBsZXhfYXJyYXknKS5Db21wbGV4QXJyYXk7XG5cbi8vIG1vZGlmaWVzIENvbXBsZXhBcnJheVxudmFyIGZmdCA9IHJlcXVpcmUoJy4uL2xpYi9qc2ZmdC9mZnQnKTtcblxubW9kdWxlLmV4cG9ydHMgPSBmdW5jdGlvbihhdWRpb0NvbnRleHQsc3JjLGJ1ZlNpemUsY2FsbGJhY2spe1xuXHRcblx0Ly9JIGFtIG15c2VsZlxuXHR2YXIgc2VsZiA9IHRoaXM7XG5cdFxuXHRzZWxmLmZlYXR1cmVFeHRyYWN0b3JzID0gcmVxdWlyZSgnLi9leHRyYWN0b3JzJyk7XG5cblx0Ly9kZWZhdWx0IGJ1ZmZlciBzaXplXG5cdHZhciBidWZmZXJTaXplID0gYnVmU2l6ZSA/IGJ1ZlNpemUgOiAyNTY7XG5cblx0Ly9pbml0aWFsIHNvdXJjZVxuXHR2YXIgc291cmNlID0gc3JjO1xuXG5cdC8vY2FsbGJhY2sgY29udHJvbGxlcnNcblx0dmFyIEVYVFJBQ1RJT05fU1RBUlRFRCA9IGZhbHNlO1xuXHR2YXIgX2ZlYXR1cmVzVG9FeHRyYWN0O1xuXG5cdC8vdXRpbGl0aWVzXG5cdHZhciDCtSA9IGZ1bmN0aW9uKGksIGFtcGxpdHVkZVNwZWN0KXtcblx0XHR2YXIgbnVtZXJhdG9yID0gMDtcblx0XHR2YXIgZGVub21pbmF0b3IgPSAwO1xuXHRcdGZvcih2YXIgayA9IDA7IGsgPCBhbXBsaXR1ZGVTcGVjdC5sZW5ndGg7IGsrKyl7XG5cdFx0XHRudW1lcmF0b3IgKz0gTWF0aC5wb3coayxpKSpNYXRoLmFicyhhbXBsaXR1ZGVTcGVjdFtrXSk7XG5cdFx0XHRkZW5vbWluYXRvciArPSBhbXBsaXR1ZGVTcGVjdFtrXTtcblx0XHR9XG5cdFx0cmV0dXJuIG51bWVyYXRvci9kZW5vbWluYXRvcjtcblx0fTtcblxuXHR2YXIgaXNQb3dlck9mVHdvID0gZnVuY3Rpb24obnVtKSB7XG5cdFx0d2hpbGUgKCgobnVtICUgMikgPT0gMCkgJiYgbnVtID4gMSkge1xuXHRcdFx0bnVtIC89IDI7XG5cdFx0fVxuXHRcdHJldHVybiAobnVtID09IDEpO1xuXHR9O1xuXG5cdC8vaW5pdGlsaXplIGJhcmsgc2NhbGUgKHRvIGJlIHVzZWQgaW4gbW9zdCBwZXJjZXB0dWFsIGZlYXR1cmVzKS5cblx0c2VsZi5iYXJrU2NhbGUgPSBuZXcgRmxvYXQzMkFycmF5KGJ1ZlNpemUpO1xuXG5cdGZvcih2YXIgaSA9IDA7IGkgPCBzZWxmLmJhcmtTY2FsZS5sZW5ndGg7IGkrKyl7XG5cdFx0c2VsZi5iYXJrU2NhbGVbaV0gPSBpKmF1ZGlvQ29udGV4dC5zYW1wbGVSYXRlLyhidWZTaXplKTtcblx0XHRzZWxmLmJhcmtTY2FsZVtpXSA9IDEzKk1hdGguYXRhbihzZWxmLmJhcmtTY2FsZVtpXS8xMzE1LjgpICsgMy41KiBNYXRoLmF0YW4oTWF0aC5wb3coKHNlbGYuYmFya1NjYWxlW2ldLzc1MTgpLDIpKTtcblx0fVxuXG5cdC8vV0lORE9XSU5HXG5cdC8vc2V0IGRlZmF1bHRcblx0c2VsZi53aW5kb3dpbmdGdW5jdGlvbiA9IFwiaGFubmluZ1wiO1xuXG5cdC8vY3JlYXRlIHdpbmRvd3Ncblx0c2VsZi5oYW5uaW5nID0gbmV3IEZsb2F0MzJBcnJheShidWZTaXplKTtcblx0Zm9yICh2YXIgaSA9IDA7IGkgPCBidWZTaXplOyBpKyspIHtcblx0XHQvL0FjY29yZGluZyB0byB0aGUgUiBkb2N1bWVudGF0aW9uIGh0dHA6Ly9yZ20ub2dhbGFiLm5ldC9SR00vUl9yZGZpbGU/Zj1HRU5FQXJlYWQvbWFuL2hhbm5pbmcud2luZG93LlJkJmQ9Ul9DQ1xuXHRcdHNlbGYuaGFubmluZ1tpXSA9IDAuNSAtIDAuNSpNYXRoLmNvcygyKk1hdGguUEkqaS8oYnVmU2l6ZS0xKSk7XG5cdH1cblxuXHRzZWxmLmhhbW1pbmcgPSBuZXcgRmxvYXQzMkFycmF5KGJ1ZlNpemUpO1xuXHRmb3IgKHZhciBpID0gMDsgaSA8IGJ1ZlNpemU7IGkrKykge1xuXHRcdC8vQWNjb3JkaW5nIHRvIGh0dHA6Ly91ay5tYXRod29ya3MuY29tL2hlbHAvc2lnbmFsL3JlZi9oYW1taW5nLmh0bWxcblx0XHRzZWxmLmhhbW1pbmdbaV0gPSAwLjU0IC0gMC40NipNYXRoLmNvcygyKk1hdGguUEkqKGkvYnVmU2l6ZS0xKSk7XG5cdH1cblxuXHQvL1VORklOSVNIRUQgLSBibGFja21hbiB3aW5kb3cgaW1wbGVtZW50YXRpb25cblxuXHQvKnNlbGYuYmxhY2ttYW4gPSBuZXcgRmxvYXQzMkFycmF5KGJ1ZlNpemUpO1xuXHQvL0FjY29yZGluZyB0byBodHRwOi8vdWsubWF0aHdvcmtzLmNvbS9oZWxwL3NpZ25hbC9yZWYvYmxhY2ttYW4uaHRtbFxuXHQvL2ZpcnN0IGhhbGYgb2YgdGhlIHdpbmRvd1xuXHRmb3IgKHZhciBpID0gMDsgaSA8IChidWZTaXplICUgMikgPyAoYnVmU2l6ZSsxKS8yIDogYnVmU2l6ZS8yOyBpKyspIHtcblx0XHRzZWxmLmJsYWNrbWFuW2ldID0gMC40MiAtIDAuNSpNYXRoLmNvcygyKk1hdGguUEkqaS8oYnVmU2l6ZS0xKSkgKyAwLjA4Kk1hdGguY29zKDQqTWF0aC5QSSppLyhidWZTaXplLTEpKTtcblx0fVxuXHQvL3NlY29uZCBoYWxmIG9mIHRoZSB3aW5kb3dcblx0Zm9yICh2YXIgaSA9IGJ1ZlNpemUvMjsgaSA+IDA7IGktLSkge1xuXHRcdHNlbGYuYmxhY2ttYW5bYnVmU2l6ZSAtIGldID0gc2VsZi5ibGFja21hbltpXTtcblx0fSovXG5cblx0c2VsZi53aW5kb3dpbmcgPSBmdW5jdGlvbihzaWcsIHR5cGUpe1xuXHRcdHZhciB3aW5kb3dlZCA9IG5ldyBGbG9hdDMyQXJyYXkoc2lnLmxlbmd0aCk7XG5cdFx0dmFyIGkgLGxlbiA9IHNpZy5sZW5ndGg7XG5cblx0XHRpZiAodHlwZSA9PSBcImhhbm5pbmdcIikge1xuXHRcdFx0Zm9yIChpID0gMDsgaSA8IGxlbjsgaSsrKSB7XG5cdFx0XHRcdHdpbmRvd2VkW2ldID0gc2lnW2ldKnNlbGYuaGFubmluZ1tpXTtcblx0XHRcdH1cblx0XHR9XG5cdFx0ZWxzZSBpZiAodHlwZSA9PSBcImhhbW1pbmdcIikge1xuXHRcdFx0Zm9yIChpID0gMDsgaSA8IGxlbjsgaSsrKSB7XG5cdFx0XHRcdHdpbmRvd2VkW2ldID0gc2lnW2ldKnNlbGYuaGFtbWluZ1tpXTtcblx0XHRcdH1cblx0XHR9XG5cdFx0ZWxzZSBpZiAodHlwZSA9PSBcImJsYWNrbWFuXCIpIHtcblx0XHRcdGZvciAoaSA9IDA7IGkgPCBsZW47IGkrKykge1xuXHRcdFx0XHR3aW5kb3dlZFtpXSA9IHNpZ1tpXSpzZWxmLmJsYWNrbWFuW2ldO1xuXHRcdFx0fVxuXHRcdH1cblxuXHRcdHJldHVybiB3aW5kb3dlZDtcblx0fTtcblxuXHQvL3NvdXJjZSBzZXR0ZXIgbWV0aG9kXG5cdHNlbGYuc2V0U291cmNlID0gZnVuY3Rpb24oX3NyYykge1xuXHRcdHNvdXJjZSA9IF9zcmM7XG5cdFx0c291cmNlLmNvbm5lY3Qod2luZG93LnNwbik7XG5cdH07XG5cblxuXG5cdGlmIChpc1Bvd2VyT2ZUd28oYnVmZmVyU2l6ZSkgJiYgYXVkaW9Db250ZXh0KSB7XG5cdFx0XHRzZWxmLmZlYXR1cmVJbmZvID0ge1xuXHRcdFx0XHRcImJ1ZmZlclwiOiB7XG5cdFx0XHRcdFx0XCJ0eXBlXCI6IFwiYXJyYXlcIlxuXHRcdFx0XHR9LFxuXHRcdFx0XHRcInJtc1wiOiB7XG5cdFx0XHRcdFx0XCJ0eXBlXCI6IFwibnVtYmVyXCJcblx0XHRcdFx0fSxcblx0XHRcdFx0XCJlbmVyZ3lcIjoge1xuXHRcdFx0XHRcdFwidHlwZVwiOiBcIm51bWJlclwiXG5cdFx0XHRcdH0sXG5cdFx0XHRcdFwiemNyXCI6IHtcblx0XHRcdFx0XHRcInR5cGVcIjogXCJudW1iZXJcIlxuXHRcdFx0XHR9LFxuXHRcdFx0XHRcImNvbXBsZXhTcGVjdHJ1bVwiOiB7XG5cdFx0XHRcdFx0XCJ0eXBlXCI6IFwibXVsdGlwbGVBcnJheXNcIixcblx0XHRcdFx0XHRcImFycmF5TmFtZXNcIjoge1xuXHRcdFx0XHRcdFx0XCIxXCI6IFwicmVhbFwiLFxuXHRcdFx0XHRcdFx0XCIyXCI6IFwiaW1hZ1wiXG5cdFx0XHRcdFx0fVxuXHRcdFx0XHR9LFxuXHRcdFx0XHRcImFtcGxpdHVkZVNwZWN0cnVtXCI6IHtcblx0XHRcdFx0XHRcInR5cGVcIjogXCJhcnJheVwiXG5cdFx0XHRcdH0sXG5cdFx0XHRcdFwicG93ZXJTcGVjdHJ1bVwiOiB7XG5cdFx0XHRcdFx0XCJ0eXBlXCI6IFwiYXJyYXlcIlxuXHRcdFx0XHR9LFxuXHRcdFx0XHRcInNwZWN0cmFsQ2VudHJvaWRcIjoge1xuXHRcdFx0XHRcdFwidHlwZVwiOiBcIm51bWJlclwiXG5cdFx0XHRcdH0sXG5cdFx0XHRcdFwic3BlY3RyYWxGbGF0bmVzc1wiOiB7XG5cdFx0XHRcdFx0XCJ0eXBlXCI6IFwibnVtYmVyXCJcblx0XHRcdFx0fSxcblx0XHRcdFx0XCJzcGVjdHJhbFNsb3BlXCI6IHtcblx0XHRcdFx0XHRcInR5cGVcIjogXCJudW1iZXJcIlxuXHRcdFx0XHR9LFxuXHRcdFx0XHRcInNwZWN0cmFsUm9sbG9mZlwiOiB7XG5cdFx0XHRcdFx0XCJ0eXBlXCI6IFwibnVtYmVyXCJcblx0XHRcdFx0fSxcblx0XHRcdFx0XCJzcGVjdHJhbFNwcmVhZFwiOiB7XG5cdFx0XHRcdFx0XCJ0eXBlXCI6IFwibnVtYmVyXCJcblx0XHRcdFx0fSxcblx0XHRcdFx0XCJzcGVjdHJhbFNrZXduZXNzXCI6IHtcblx0XHRcdFx0XHRcInR5cGVcIjogXCJudW1iZXJcIlxuXHRcdFx0XHR9LFxuXHRcdFx0XHRcInNwZWN0cmFsS3VydG9zaXNcIjoge1xuXHRcdFx0XHRcdFwidHlwZVwiOiBcIm51bWJlclwiXG5cdFx0XHRcdH0sXG5cdFx0XHRcdFwibG91ZG5lc3NcIjoge1xuXHRcdFx0XHRcdFwidHlwZVwiOiBcIm11bHRpcGxlQXJyYXlzXCIsXG5cdFx0XHRcdFx0XCJhcnJheU5hbWVzXCI6IHtcblx0XHRcdFx0XHRcdFwiMVwiOiBcInRvdGFsXCIsXG5cdFx0XHRcdFx0XHRcIjJcIjogXCJzcGVjaWZpY1wiXG5cdFx0XHRcdFx0fVxuXHRcdFx0XHR9LFxuXHRcdFx0XHRcInBlcmNlcHR1YWxTcHJlYWRcIjoge1xuXHRcdFx0XHRcdFwidHlwZVwiOiBcIm51bWJlclwiXG5cdFx0XHRcdH0sXG5cdFx0XHRcdFwicGVyY2VwdHVhbFNoYXJwbmVzc1wiOiB7XG5cdFx0XHRcdFx0XCJ0eXBlXCI6IFwibnVtYmVyXCJcblx0XHRcdFx0fSxcblx0XHRcdFx0XCJtZmNjXCI6IHtcblx0XHRcdFx0XHRcInR5cGVcIjogXCJhcnJheVwiXG5cdFx0XHRcdH1cblx0XHRcdH07XG5cblx0XHRcdC8vY3JlYXRlIGNvbXBsZXhhcnJheSB0byBob2xkIHRoZSBzcGVjdHJ1bVxuXHRcdFx0dmFyIGRhdGEgPSBuZXcgQ29tcGxleEFycmF5KGJ1ZmZlclNpemUpO1xuXG5cdFx0XHQvL2NyZWF0ZSBub2Rlc1xuXHRcdFx0d2luZG93LnNwbiA9IGF1ZGlvQ29udGV4dC5jcmVhdGVTY3JpcHRQcm9jZXNzb3IoYnVmZmVyU2l6ZSwxLDEpO1xuXHRcdFx0c3BuLmNvbm5lY3QoYXVkaW9Db250ZXh0LmRlc3RpbmF0aW9uKTtcblxuXHRcdFx0d2luZG93LnNwbi5vbmF1ZGlvcHJvY2VzcyA9IGZ1bmN0aW9uKGUpIHtcblx0XHRcdFx0Ly90aGlzIGlzIHRvIG9idGFpbiB0aGUgY3VycmVudCBhbXBsaXR1ZGUgc3BlY3RydW1cblx0XHRcdFx0dmFyIGlucHV0RGF0YSA9IGUuaW5wdXRCdWZmZXIuZ2V0Q2hhbm5lbERhdGEoMCk7XG5cdFx0XHRcdHNlbGYuc2lnbmFsID0gaW5wdXREYXRhO1xuXHRcdFx0XHR2YXIgd2luZG93ZWRTaWduYWwgPSBzZWxmLndpbmRvd2luZyhzZWxmLnNpZ25hbCwgc2VsZi53aW5kb3dpbmdGdW5jdGlvbik7XG5cblx0XHRcdFx0Ly9tYXAgdGltZSBkb21haW5cblx0XHRcdFx0ZGF0YS5tYXAoZnVuY3Rpb24odmFsdWUsIGksIG4pIHtcblx0XHRcdFx0XHR2YWx1ZS5yZWFsID0gd2luZG93ZWRTaWduYWxbaV07XG5cdFx0XHRcdH0pO1xuXHRcdFx0XHQvL3RyYW5zZm9ybVxuXHRcdFx0XHR2YXIgc3BlYyA9IGRhdGEuRkZUKCk7XG5cdFx0XHRcdC8vYXNzaWduIHRvIG1leWRhXG5cdFx0XHRcdHNlbGYuY29tcGxleFNwZWN0cnVtID0gc3BlYztcblx0XHRcdFx0c2VsZi5hbXBTcGVjdHJ1bSA9IG5ldyBGbG9hdDMyQXJyYXkoYnVmZmVyU2l6ZS8yKTtcblx0XHRcdFx0Ly9jYWxjdWxhdGUgYW1wbGl0dWRlXG5cdFx0XHRcdGZvciAodmFyIGkgPSAwOyBpIDwgYnVmZmVyU2l6ZS8yOyBpKyspIHtcblx0XHRcdFx0XHRzZWxmLmFtcFNwZWN0cnVtW2ldID0gTWF0aC5zcXJ0KE1hdGgucG93KHNwZWMucmVhbFtpXSwyKSArIE1hdGgucG93KHNwZWMuaW1hZ1tpXSwyKSk7XG5cblx0XHRcdFx0fVxuXHRcdFx0XHQvL2NhbGwgY2FsbGJhY2sgaWYgYXBwbGljYWJsZVxuXHRcdFx0XHRpZiAodHlwZW9mIGNhbGxiYWNrID09PSBcImZ1bmN0aW9uXCIgJiYgRVhUUkFDVElPTl9TVEFSVEVEKSB7XG5cdFx0XHRcdFx0Y2FsbGJhY2soc2VsZi5nZXQoX2ZlYXR1cmVzVG9FeHRyYWN0KSk7XG5cdFx0XHRcdH1cblxuXHRcdFx0fTtcblxuXHRcdFx0c2VsZi5zdGFydCA9IGZ1bmN0aW9uKGZlYXR1cmVzKSB7XG5cdFx0XHRcdF9mZWF0dXJlc1RvRXh0cmFjdCA9IGZlYXR1cmVzO1xuXHRcdFx0XHRFWFRSQUNUSU9OX1NUQVJURUQgPSB0cnVlO1xuXHRcdFx0fTtcblxuXHRcdFx0c2VsZi5zdG9wID0gZnVuY3Rpb24oKSB7XG5cdFx0XHRcdEVYVFJBQ1RJT05fU1RBUlRFRCA9IGZhbHNlO1xuXHRcdFx0fTtcblxuXHRcdFx0c2VsZi5hdWRpb0NvbnRleHQgPSBhdWRpb0NvbnRleHQ7XG5cblx0XHRcdHNlbGYuZ2V0ID0gZnVuY3Rpb24oZmVhdHVyZSkge1xuXHRcdFx0XHRpZih0eXBlb2YgZmVhdHVyZSA9PT0gXCJvYmplY3RcIil7XG5cdFx0XHRcdFx0dmFyIHJlc3VsdHMgPSB7fTtcblx0XHRcdFx0XHRmb3IgKHZhciB4ID0gMDsgeCA8IGZlYXR1cmUubGVuZ3RoOyB4Kyspe1xuXHRcdFx0XHRcdFx0dHJ5e1xuXHRcdFx0XHRcdFx0XHRyZXN1bHRzW2ZlYXR1cmVbeF1dID0gKHNlbGYuZmVhdHVyZUV4dHJhY3RvcnNbZmVhdHVyZVt4XV0oYnVmZmVyU2l6ZSwgc2VsZikpO1xuXHRcdFx0XHRcdFx0fSBjYXRjaCAoZSl7XG5cdFx0XHRcdFx0XHRcdGNvbnNvbGUuZXJyb3IoZSk7XG5cdFx0XHRcdFx0XHR9XG5cdFx0XHRcdFx0fVxuXHRcdFx0XHRcdHJldHVybiByZXN1bHRzO1xuXHRcdFx0XHR9XG5cdFx0XHRcdGVsc2UgaWYgKHR5cGVvZiBmZWF0dXJlID09PSBcInN0cmluZ1wiKXtcblx0XHRcdFx0XHRyZXR1cm4gc2VsZi5mZWF0dXJlRXh0cmFjdG9yc1tmZWF0dXJlXShidWZmZXJTaXplLCBzZWxmKTtcblx0XHRcdFx0fVxuXHRcdFx0XHRlbHNle1xuXHRcdFx0XHRcdHRocm93IFwiSW52YWxpZCBGZWF0dXJlIEZvcm1hdFwiO1xuXHRcdFx0XHR9XG5cdFx0XHR9O1xuXHRcdFx0c291cmNlLmNvbm5lY3Qod2luZG93LnNwbiwgMCwgMCk7XG5cdFx0XHRyZXR1cm4gc2VsZjtcblx0fVxuXHRlbHNlIHtcblx0XHQvL2hhbmRsZSBlcnJvcnNcblx0XHRpZiAodHlwZW9mIGF1ZGlvQ29udGV4dCA9PSBcInVuZGVmaW5lZFwiKSB7XG5cdFx0XHR0aHJvdyBcIkF1ZGlvQ29udGV4dCB3YXNuJ3Qgc3BlY2lmaWVkOiBNZXlkYSB3aWxsIG5vdCBydW4uXCI7XG5cdFx0fVxuXHRcdGVsc2Uge1xuXHRcdFx0dGhyb3cgXCJCdWZmZXIgc2l6ZSBpcyBub3QgYSBwb3dlciBvZiB0d286IE1leWRhIHdpbGwgbm90IHJ1bi5cIjtcblx0XHR9XG5cdH1cbn07Il19

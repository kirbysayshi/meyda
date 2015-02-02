!function(e){if("object"==typeof exports&&"undefined"!=typeof module)module.exports=e();else if("function"==typeof define&&define.amd)define([],e);else{var f;"undefined"!=typeof window?f=window:"undefined"!=typeof global?f=global:"undefined"!=typeof self&&(f=self),f.Meyda=e()}}(function(){var define,module,exports;return (function e(t,n,r){function s(o,u){if(!n[o]){if(!t[o]){var a=typeof require=="function"&&require;if(!u&&a)return a(o,!0);if(i)return i(o,!0);var f=new Error("Cannot find module '"+o+"'");throw f.code="MODULE_NOT_FOUND",f}var l=n[o]={exports:{}};t[o][0].call(l.exports,function(e){var n=t[o][1][e];return s(n?n:e)},l,l.exports,e,t,n,r)}return n[o].exports}var i=typeof require=="function"&&require;for(var o=0;o<r.length;o++)s(r[o]);return s})({1:[function(require,module,exports){
module.exports = function() {
  "use strict";
  'use strict';
  !function(exports, undefined) {
    var DefaultArrayType = Float32Array,
        sqrt = Math.sqrt,
        sqr = function(number) {
          return Math.pow(number, 2);
        },
        isComplexArray,
        ComplexArray;
    exports.isComplexArray = isComplexArray = function(obj) {
      return obj !== undefined && obj.hasOwnProperty !== undefined && obj.hasOwnProperty('real') && obj.hasOwnProperty('imag');
    };
    exports.ComplexArray = ComplexArray = function(other, opt_array_type) {
      if (isComplexArray(other)) {
        this.ArrayType = other.ArrayType;
        this.real = new this.ArrayType(other.real);
        this.imag = new this.ArrayType(other.imag);
      } else {
        this.ArrayType = opt_array_type || DefaultArrayType;
        this.real = new this.ArrayType(other);
        this.imag = new this.ArrayType(this.real.length);
      }
      this.length = this.real.length;
    };
    ComplexArray.prototype.toString = function() {
      var components = [];
      this.forEach(function(c_value, i) {
        components.push('(' + c_value.real.toFixed(2) + ',' + c_value.imag.toFixed(2) + ')');
      });
      return '[' + components.join(',') + ']';
    };
    ComplexArray.prototype.map = function(mapper) {
      var i,
          n = this.length,
          c_value = {};
      for (i = 0; i < n; i++) {
        c_value.real = this.real[i];
        c_value.imag = this.imag[i];
        mapper(c_value, i, n);
        this.real[i] = c_value.real;
        this.imag[i] = c_value.imag;
      }
      return this;
    };
    ComplexArray.prototype.forEach = function(iterator) {
      var i,
          n = this.length,
          c_value = {};
      for (i = 0; i < n; i++) {
        c_value.real = this.real[i];
        c_value.imag = this.imag[i];
        iterator(c_value, i, n);
      }
    };
    ComplexArray.prototype.conjugate = function() {
      return (new ComplexArray(this)).map(function(value) {
        value.imag *= -1;
      });
    };
    function iterable(obj) {
      if (!obj.forEach)
        obj.forEach = function(iterator) {
          var i,
              n = this.length;
          for (i = 0; i < n; i++)
            iterator(this[i], i, n);
        };
      return obj;
    }
    ComplexArray.prototype.magnitude = function() {
      var mags = new this.ArrayType(this.length);
      this.forEach(function(value, i) {
        mags[i] = sqrt(sqr(value.real) + sqr(value.imag));
      });
      return iterable(mags);
    };
  }(typeof exports === 'undefined' && (this.complex_array = {}) || exports);
  return {};
}.call(Reflect.global);


//# sourceURL=/Volumes/Home/Dev/meyda-modular/lib/jsfft/complex_array.js
},{}],2:[function(require,module,exports){
module.exports = function() {
  "use strict";
  'use strict';
  !function(exports, complex_array) {
    var ComplexArray = complex_array.ComplexArray,
        PI = Math.PI,
        SQRT1_2 = Math.SQRT1_2,
        sqrt = Math.sqrt,
        cos = Math.cos,
        sin = Math.sin;
    ComplexArray.prototype.FFT = function() {
      return FFT(this, false);
    };
    exports.FFT = function(input) {
      return ensureComplexArray(input).FFT();
    };
    ComplexArray.prototype.InvFFT = function() {
      return FFT(this, true);
    };
    exports.InvFFT = function(input) {
      return ensureComplexArray(input).InvFFT();
    };
    ComplexArray.prototype.frequencyMap = function(filterer) {
      return this.FFT().map(filterer).InvFFT();
    };
    exports.frequencyMap = function(input, filterer) {
      return ensureComplexArray(input).frequencyMap(filterer);
    };
    function ensureComplexArray(input) {
      return complex_array.isComplexArray(input) && input || new ComplexArray(input);
    }
    function FFT(input, inverse) {
      var n = input.length;
      if (n & (n - 1)) {
        return FFT_Recursive(input, inverse);
      } else {
        return FFT_2_Iterative(input, inverse);
      }
    }
    function FFT_Recursive(input, inverse) {
      var n = input.length,
          i,
          j,
          output,
          f_r,
          f_i,
          del_f_r,
          del_f_i,
          p,
          m,
          normalisation,
          recursive_result,
          _swap,
          _real,
          _imag;
      if (n === 1) {
        return input;
      }
      output = new ComplexArray(n, input.ArrayType);
      p = LowestOddFactor(n);
      m = n / p;
      normalisation = 1 / sqrt(p);
      recursive_result = new ComplexArray(m, input.ArrayType);
      for (j = 0; j < p; j++) {
        for (i = 0; i < m; i++) {
          recursive_result.real[i] = input.real[i * p + j];
          recursive_result.imag[i] = input.imag[i * p + j];
        }
        if (m > 1) {
          recursive_result = FFT(recursive_result, inverse);
        }
        del_f_r = cos(2 * PI * j / n);
        del_f_i = (inverse ? -1 : 1) * sin(2 * PI * j / n);
        f_r = 1;
        f_i = 0;
        for (i = 0; i < n; i++) {
          _real = recursive_result.real[i % m];
          _imag = recursive_result.imag[i % m];
          output.real[i] += f_r * _real - f_i * _imag;
          output.imag[i] += f_r * _imag + f_i * _real;
          _swap = f_r * del_f_r - f_i * del_f_i;
          f_i = f_r * del_f_i + f_i * del_f_r;
          f_r = _swap;
        }
      }
      for (i = 0; i < n; i++) {
        input.real[i] = normalisation * output.real[i];
        input.imag[i] = normalisation * output.imag[i];
      }
      return input;
    }
    function FFT_2_Iterative(input, inverse) {
      var n = input.length,
          i,
          j,
          output,
          output_r,
          output_i,
          f_r,
          f_i,
          del_f_r,
          del_f_i,
          temp,
          l_index,
          r_index,
          left_r,
          left_i,
          right_r,
          right_i,
          width;
      output = BitReverseComplexArray(input);
      output_r = output.real;
      output_i = output.imag;
      width = 1;
      while (width < n) {
        del_f_r = cos(PI / width);
        del_f_i = (inverse ? -1 : 1) * sin(PI / width);
        for (i = 0; i < n / (2 * width); i++) {
          f_r = 1;
          f_i = 0;
          for (j = 0; j < width; j++) {
            l_index = 2 * i * width + j;
            r_index = l_index + width;
            left_r = output_r[l_index];
            left_i = output_i[l_index];
            right_r = f_r * output_r[r_index] - f_i * output_i[r_index];
            right_i = f_i * output_r[r_index] + f_r * output_i[r_index];
            output_r[l_index] = SQRT1_2 * (left_r + right_r);
            output_i[l_index] = SQRT1_2 * (left_i + right_i);
            output_r[r_index] = SQRT1_2 * (left_r - right_r);
            output_i[r_index] = SQRT1_2 * (left_i - right_i);
            temp = f_r * del_f_r - f_i * del_f_i;
            f_i = f_r * del_f_i + f_i * del_f_r;
            f_r = temp;
          }
        }
        width <<= 1;
      }
      return output;
    }
    function BitReverseIndex(index, n) {
      var bitreversed_index = 0;
      while (n > 1) {
        bitreversed_index <<= 1;
        bitreversed_index += index & 1;
        index >>= 1;
        n >>= 1;
      }
      return bitreversed_index;
    }
    function BitReverseComplexArray(array) {
      var n = array.length,
          flips = {},
          swap,
          i;
      for (i = 0; i < n; i++) {
        var r_i = BitReverseIndex(i, n);
        if (flips.hasOwnProperty(i) || flips.hasOwnProperty(r_i))
          continue;
        swap = array.real[r_i];
        array.real[r_i] = array.real[i];
        array.real[i] = swap;
        swap = array.imag[r_i];
        array.imag[r_i] = array.imag[i];
        array.imag[i] = swap;
        flips[i] = flips[r_i] = true;
      }
      return array;
    }
    function LowestOddFactor(n) {
      var factor = 3,
          sqrt_n = sqrt(n);
      while (factor <= sqrt_n) {
        if (n % factor === 0)
          return factor;
        factor = factor + 2;
      }
      return n;
    }
  }(typeof exports === 'undefined' && (this.fft = {}) || exports, typeof require === 'undefined' && (this.complex_array) || require('./complex_array'));
  return {};
}.call(Reflect.global);


//# sourceURL=/Volumes/Home/Dev/meyda-modular/lib/jsfft/fft.js
},{"./complex_array":1}],3:[function(require,module,exports){
"use strict";
module.exports = {loudness: require('./loudness')};


//# sourceURL=/Volumes/Home/Dev/meyda-modular/src/extractors/index.js
},{"./loudness":4}],4:[function(require,module,exports){
"use strict";
function Loudness(opts) {
  if (!(this instanceof Loudness))
    return new Loudness(opts);
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
  this.process = this.process.bind(this);
  this.barkScale = opts.barkScale;
  this.bbLimits = this.computeBarkBandLimits(this.barkScale, this.normalisedSpectrum.length, this.NUM_BARK_BANDS);
}
Loudness.prototype.computeBarkBandLimits = function(barkScale, nSpectrumLength, num_bark_bands) {
  barkScale = barkScale || this.barkScale;
  nSpectrumLength = nSpectrumLength || this.normalisedSpectrum.length;
  num_bark_bands = num_bark_bands || this.NUM_BARK_BANDS;
  var currentBandEnd = barkScale[nSpectrumLength - 1] / num_bark_bands;
  var currentBand = 1;
  var bbLimits = new Int32Array(num_bark_bands + 1);
  bbLimits[0] = 0;
  for (var i = 0; i < nSpectrumLength; i++) {
    while (barkScale[i] > currentBandEnd) {
      bbLimits[currentBand++] = i;
      currentBandEnd = currentBand * barkScale[nSpectrumLength - 1] / num_bark_bands;
    }
  }
  bbLimits[num_bark_bands] = nSpectrumLength - 1;
  return bbLimits;
};
Loudness.prototype.computeSpecificLoudness = function(bbLimits, specific, nSpectrum, num_bark_bands) {
  bbLimits = bbLimits || this.bbLimits;
  specific = specific || this.specific;
  nSpectrum = nSpectrum || this.normalisedSpectrum;
  num_bark_bands = num_bark_bands || this.NUM_BARK_BANDS;
  for (var i = 0; i < num_bark_bands; i++) {
    var sum = 0;
    for (var j = bbLimits[i]; j < bbLimits[i + 1]; j++) {
      sum += nSpectrum[j];
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
Loudness.prototype.process = function() {
  var barkScale = this.barkScale;
  var nSpectrum = this.normalisedSpectrum;
  var nSpectrumLength = nSpectrum.length;
  var NUM_BARK_BANDS = this.NUM_BARK_BANDS;
  var specific = this.specific;
  var bbLimits = this.bbLimits;
  var spec = this.computeSpecificLoudness(bbLimits, specific, nSpectrum, NUM_BARK_BANDS);
  var total = this.sumArray(spec);
  return {
    specific: spec,
    total: total
  };
};
module.exports = Loudness;


//# sourceURL=/Volumes/Home/Dev/meyda-modular/src/extractors/loudness.js
},{}],5:[function(require,module,exports){
"use strict";
module.exportts = {
  "buffer": {"type": "array"},
  "rms": {"type": "number"},
  "energy": {"type": "number"},
  "zcr": {"type": "number"},
  "complexSpectrum": {
    "type": "multipleArrays",
    "arrayNames": {
      "1": "real",
      "2": "imag"
    }
  },
  "amplitudeSpectrum": {"type": "array"},
  "powerSpectrum": {"type": "array"},
  "spectralCentroid": {"type": "number"},
  "spectralFlatness": {"type": "number"},
  "spectralSlope": {"type": "number"},
  "spectralRolloff": {"type": "number"},
  "spectralSpread": {"type": "number"},
  "spectralSkewness": {"type": "number"},
  "spectralKurtosis": {"type": "number"},
  "loudness": {
    "type": "multipleArrays",
    "arrayNames": {
      "1": "total",
      "2": "specific"
    }
  },
  "perceptualSpread": {"type": "number"},
  "perceptualSharpness": {"type": "number"},
  "mfcc": {"type": "array"}
};


//# sourceURL=/Volumes/Home/Dev/meyda-modular/src/feature-info.js
},{}],6:[function(require,module,exports){
"use strict";
var ComplexArray = require('../lib/jsfft/complex_array').ComplexArray,
    fft = require('../lib/jsfft/fft');
var $__3 = require('./utils'),
    isPowerOfTwo = $__3.isPowerOfTwo,
    µ = $__3.µ;
var Meyda = function Meyda(audioContext, src, bufSize, callback) {
  var $__0 = this;
  if (!isPowerOfTwo(bufSize)) {
    throw new Error("Buffer size is not a power of two: Meyda will not run.");
  }
  if (!audioContext) {
    throw new Error("AudioContext wasn't specified: Meyda will not run.");
  }
  var bufferSize = bufSize || 256;
  var sampleRate = audioContext.sampleRate;
  var source = src;
  this.audioContext = audioContext;
  this.featureExtractors = {};
  this.EXTRACTION_STARTED = false;
  this._featuresToExtract = null;
  this.windowingFunction = "hanning";
  this.barkScale = this.computeBarkScale(bufferSize, sampleRate);
  this.hanning = this.computeHanning(bufferSize);
  this.hamming = this.computeHamming(bufferSize);
  this.featureInfo = require('./feature-info');
  var computedSpectrumData = this.computeSpectrumData(bufferSize);
  this.spectrumData = computedSpectrumData.data;
  this.complexSpectrum = computedSpectrumData.spectrum;
  this.ampSpectrum = computedSpectrumData.ampSpectrum;
  this.initialiseExtractors();
  window.spn = audioContext.createScriptProcessor(bufferSize, 1, 1);
  window.spn.onaudioprocess = (function(e) {
    var signal = $__0.signal = e.inputBuffer.getChannelData(0);
    var data = $__0.spectrumData;
    var spec = $__0.complexSpectrum;
    var ampSpectrum = $__0.ampSpectrum;
    var windowedSignal = $__0.computeWindow(signal, $__0.windowingFunction);
    data.map(function(value, i, n) {
      value.real = windowedSignal[i];
    });
    $__0.computeAmplitude(spec, ampSpectrum, bufferSize);
    if (typeof callback === "function" && EXTRACTION_STARTED) {
      callback($__0.get($__0._featuresToExtract));
    }
  });
  window.spn.connect(audioContext.destination);
  source.connect(window.spn, 0, 0);
};
($traceurRuntime.createClass)(Meyda, {
  computeAmplitude: function(complexSpectrum, ampSpectrum, bufferSize) {
    complexSpectrum = complexSpectrum || this.complexSpectrum;
    ampSpectrum = ampSpectrum || this.ampSpectrum;
    bufferSize = bufferSize || this.bufferSize;
    for (var i = 0; i < bufferSize / 2; i++) {
      ampSpectrum[i] = Math.sqrt(Math.pow(complexSpectrum.real[i], 2) + Math.pow(complexSpectrum.imag[i], 2));
    }
  },
  computeHamming: function(bufferSize) {
    bufferSize = bufferSize || this.bufferSize;
    var hamming = new Float32Array(bufferSize);
    for (var i = 0; i < bufferSize; i++) {
      hamming[i] = 0.54 - 0.46 * Math.cos(2 * Math.PI * (i / bufferSize - 1));
    }
    return hamming;
  },
  computeHanning: function(bufferSize) {
    bufferSize = bufferSize || this.bufferSize;
    var hanning = new Float32Array(bufferSize);
    for (var i = 0; i < bufferSize; i++) {
      hanning[i] = 0.5 - 0.5 * Math.cos(2 * Math.PI * i / (bufferSize - 1));
    }
    return hanning;
  },
  computeWindow: function(sig, type) {
    var i,
        len = sig.length;
    var windowed = new Float32Array(len);
    for (i = 0; i < len; i++) {
      windowed[i] = sig[i] * this[type][i];
    }
    return windowed;
  },
  computeBarkScale: function(bufferSize, sampleRate) {
    bufferSize = bufferSize || this.bufferSize;
    sampleRate = sampleRate || this.sampleRate;
    var barkScale = new Float32Array(bufferSize);
    for (var i = 0; i < bufferSize; i++) {
      barkScale[i] = i * sampleRate / (bufferSize);
      barkScale[i] = 13 * Math.atan(barkScale[i] / 1315.8) + 3.5 * Math.atan(Math.pow((barkScale[i] / 7518), 2));
    }
    return barkScale;
  },
  computeSpectrumData: function(bufferSize) {
    bufferSize = bufferSize || this.bufferSize;
    var data = new ComplexArray(bufferSize);
    var spectrum = data.FFT();
    var ampSpectrum = new Float32Array(bufferSize / 2);
    return {
      data: data,
      spectrum: spectrum,
      ampSpectrum: ampSpectrum
    };
  },
  initialiseExtractors: function() {
    var extractors = require('./extractors');
    var loudness = extractors.loudness({
      NUM_BARK_BANDS: 24,
      barkScale: this.barkScale,
      normalisedSpectrum: this.ampSpectrum,
      sampleRate: this.audioContext.sampleRate
    });
    this.featureExtractors.loudness = loudness;
    this.featureInfo.loudness = loudness.info;
  },
  setSource: function(_src) {
    _src.connect(window.spn);
  },
  start: function(features) {
    this._featuresToExtract = features;
    this.EXTRACTION_STARTED = true;
  },
  stop: function() {
    this._featuresToExtract = null;
    this.EXTRACTION_STARTED = false;
  },
  get: function(feature) {
    if (typeof feature === "object") {
      var results = {};
      for (var x = 0; x < feature.length; x++) {
        try {
          results[feature[x]] = (this.featureExtractors[feature[x]].process(this.signal));
        } catch (e) {
          console.error(e);
        }
      }
      return results;
    } else if (typeof feature === "string") {
      return this.featureExtractors[feature].process(this.signal);
    } else {
      throw new Error("Invalid Feature Format");
    }
  }
}, {});
module.exports = Meyda;


//# sourceURL=/Volumes/Home/Dev/meyda-modular/src/meyda.js
},{"../lib/jsfft/complex_array":1,"../lib/jsfft/fft":2,"./extractors":3,"./feature-info":5,"./utils":7}],7:[function(require,module,exports){
"use strict";
module.exports.µ = function(i, amplitudeSpect) {
  var numerator = 0;
  var denominator = 0;
  for (var k = 0; k < amplitudeSpect.length; k++) {
    numerator += Math.pow(k, i) * Math.abs(amplitudeSpect[k]);
    denominator += amplitudeSpect[k];
  }
  return numerator / denominator;
};
module.exports.isPowerOfTwo = function(num) {
  while (((num % 2) === 0) && num > 1) {
    num /= 2;
  }
  return (num == 1);
};


//# sourceURL=/Volumes/Home/Dev/meyda-modular/src/utils.js
},{}]},{},[6])(6)
});
//# sourceMappingURL=data:application/json;base64,eyJ2ZXJzaW9uIjozLCJzb3VyY2VzIjpbIi4uLy4uLy4uLy4uL3Vzci9sb2NhbC9saWIvbm9kZV9tb2R1bGVzL2Jyb3dzZXJpZnkvbm9kZV9tb2R1bGVzL2Jyb3dzZXItcGFjay9fcHJlbHVkZS5qcyIsIkB0cmFjZXVyL2dlbmVyYXRlZC9UZW1wbGF0ZVBhcnNlci80IiwiL1ZvbHVtZXMvSG9tZS9EZXYvbWV5ZGEtbW9kdWxhci9zcmMvZXh0cmFjdG9ycy9pbmRleC5qcyIsIi9Wb2x1bWVzL0hvbWUvRGV2L21leWRhLW1vZHVsYXIvc3JjL2V4dHJhY3RvcnMvbG91ZG5lc3MuanMiLCIvVm9sdW1lcy9Ib21lL0Rldi9tZXlkYS1tb2R1bGFyL3NyYy9mZWF0dXJlLWluZm8uanMiLCIvVm9sdW1lcy9Ib21lL0Rldi9tZXlkYS1tb2R1bGFyL3NyYy9tZXlkYS5qcyIsIi9Wb2x1bWVzL0hvbWUvRGV2L21leWRhLW1vZHVsYXIvc3JjL3V0aWxzLmpzIl0sIm5hbWVzIjpbXSwibWFwcGluZ3MiOiJBQUFBO0FDQUEsS0FBSyxRQUFRLEVBQUksQ0FBQSxTQUFRLEFBQUM7O0FBQTFCLGFBQVcsQ0FBQztBQUVaLEVBQUMsU0FBUyxPQUFNLENBQUcsQ0FBQSxTQUFRLENBQUc7QUFFNUIsQUFFRSxNQUFBLENBQUEsZ0JBQWUsRUFBSSxhQUFXO0FBRTlCLFdBQUcsRUFBSSxDQUFBLElBQUcsS0FBSztBQUNmLFVBQUUsRUFBSSxVQUFTLE1BQUssQ0FBRztBQUFDLGVBQU8sQ0FBQSxJQUFHLElBQUksQUFBQyxDQUFDLE1BQUssQ0FBRyxFQUFBLENBQUMsQ0FBQTtRQUFDO0FBRWxELHFCQUFhO0FBQ2IsbUJBQVcsQ0FBQTtBQUViLFVBQU0sZUFBZSxFQUFJLENBQUEsY0FBYSxFQUFJLFVBQVMsR0FBRSxDQUFHO0FBQ3RELFdBQU8sQ0FBQSxHQUFFLElBQU0sVUFBUSxDQUFBLEVBQ3JCLENBQUEsR0FBRSxlQUFlLElBQU0sVUFBUSxDQUFBLEVBQy9CLENBQUEsR0FBRSxlQUFlLEFBQUMsQ0FBQyxNQUFLLENBQUMsQ0FBQSxF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QSxDQUFDLENBQUMsQ0FBQztJQUN0RztBQUFBLEVBQ0Q7QUFFQSxlQUFhLENBQWIsVUFBZSxVQUFTLENBQUc7QUFDMUIsYUFBUyxFQUFJLENBQUEsVUFBUyxHQUFLLENBQUEsSUFBRyxXQUFXLENBQUM7QUFFMUMsQUFBSSxNQUFBLENBQUEsT0FBTSxFQUFJLElBQUksYUFBVyxBQUFDLENBQUMsVUFBUyxDQUFDLENBQUM7QUFDMUMsUUFBUyxHQUFBLENBQUEsQ0FBQSxFQUFJLEVBQUEsQ0FBRyxDQUFBLENBQUEsRUFBSSxXQUFTLENBQUcsQ0FBQSxDQUFBLEVBQUUsQ0FBRztBQUVwQyxZQUFNLENBQUUsQ0FBQSxDQUFDLEVBQUksQ0FBQSxJQUFHLEVBQUksQ0FBQSxJQUFHLEVBQUUsQ0FBQSxJQUFHLElBQUksQUFBQyxDQUFDLENBQUEsRUFBRSxDQUFBLElBQUcsR0FBRyxDQUFBLENBQUUsRUFBQyxDQUFBLEVBQUUsV0FBUyxDQUFBLENBQUUsRUFBQSxDQUFDLENBQUMsQ0FBQztJQUM5RDtBQUFBLEFBRUEsU0FBTyxRQUFNLENBQUM7RUFDZDtBQUVELGVBQWEsQ0FBYixVQUFlLFVBQVMsQ0FBRztBQUMxQixhQUFTLEVBQUksQ0FBQSxVQUFTLEdBQUssQ0FBQSxJQUFHLFdBQVcsQ0FBQztBQUUxQyxBQUFJLE1BQUEsQ0FBQSxPQUFNLEVBQUksSUFBSSxhQUFXLEFBQUMsQ0FBQyxVQUFTLENBQUMsQ0FBQztBQUMxQyxRQUFTLEdBQUEsQ0FBQSxDQUFBLEVBQUksRUFBQSxDQUFHLENBQUEsQ0FBQSxFQUFJLFdBQVMsQ0FBRyxDQUFBLENBQUEsRUFBRSxDQUFHO0FBRXBDLFlBQU0sQ0FBRSxDQUFBLENBQUMsRUFBSSxDQUFBLEdBQUUsRUFBSSxDQUFBLEdBQUUsRUFBRSxDQUFBLElBQUcsSUFBSSxBQUFDLENBQUMsQ0FBQSxFQUFFLENBQUEsSUFBRyxHQUFHLENBQUEsQ0FBRSxFQUFBLENBQUEsQ0FBRSxFQUFDLFVBQVMsRUFBRSxFQUFBLENBQUMsQ0FBQyxDQUFDO0lBQzVEO0FBQUEsQUFFQSxTQUFPLFFB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QUM1RCxLQUFNO0FBQ0wsVUFBTSxJQUFJLE1BQUksQUFBQyxDQUFDLHdCQUF1QixDQUFDLENBQUM7SUFDMUM7QUFBQSxFQUNEO0FBQUEsS0FwUW9GO0FBd1FyRixLQUFLLFFBQVEsRUFBSSxNQUFJLENBQUM7QUFDdEI7Ozs7QUN6UUE7QUFBQSxLQUFLLFFBQVEsRUFBRSxFQUFJLFVBQVMsQ0FBQSxDQUFHLENBQUEsY0FBYSxDQUFFO0FBQzVDLEFBQUksSUFBQSxDQUFBLFNBQVEsRUFBSSxFQUFBLENBQUM7QUFDakIsQUFBSSxJQUFBLENBQUEsV0FBVSxFQUFJLEVBQUEsQ0FBQztBQUVuQixNQUFRLEdBQUEsQ0FBQSxDQUFBLEVBQUksRUFBQSxDQUFHLENBQUEsQ0FBQSxFQUFJLENBQUEsY0FBYSxPQUFPLENBQUcsQ0FBQSxDQUFBLEVBQUUsQ0FBRTtBQUM1QyxZQUFRLEdBQUssQ0FBQSxJQUFHLElBQUksQUFBQyxDQUFDLENBQUEsQ0FBRSxFQUFBLENBQUMsQ0FBQSxDQUFFLENBQUEsSUFBRyxJQUFJLEFBQUMsQ0FBQyxjQUFhLENBQUUsQ0FBQSxDQUFDLENBQUMsQ0FBQztBQUN0RCxjQUFVLEdBQUssQ0FBQSxjQUFhLENBQUUsQ0FBQSxDQUFDLENBQUM7RUFDbEM7QUFBQSxBQUVBLE9BQU8sQ0FBQSxTQUFRLEVBQUUsWUFBVSxDQUFDO0FBQzlCLENBQUM7QUFFRCxLQUFLLFFBQVEsYUFBYSxFQUFJLFVBQVMsR0FBRSxDQUFHO0FBQzFDLFFBQU8sQ0FBQyxDQUFDLEdBQUUsRUFBSSxFQUFBLENBQUMsSUFBTSxFQUFBLENBQUMsR0FBSyxDQUFBLEdBQUUsRUFBSSxFQUFBLENBQUc7QUFDbkMsTUFBRSxHQUFLLEVBQUEsQ0FBQztFQUNWO0FBQUEsQUFFQSxPQUFPLEVBQUMsR0FBRSxHQUFLLEVBQUEsQ0FBQyxDQUFDO0FBQ25CLENBQUM7QUFBQSIsImZpbGUiOiJnZW5lcmF0ZWQuanMiLCJzb3VyY2VSb290IjoiIiwic291cmNlc0NvbnRlbnQiOlsiKGZ1bmN0aW9uIGUodCxuLHIpe2Z1bmN0aW9uIHMobyx1KXtpZighbltvXSl7aWYoIXRbb10pe3ZhciBhPXR5cGVvZiByZXF1aXJlPT1cImZ1bmN0aW9uXCImJnJlcXVpcmU7aWYoIXUmJmEpcmV0dXJuIGEobywhMCk7aWYoaSlyZXR1cm4gaShvLCEwKTt2YXIgZj1uZXcgRXJyb3IoXCJDYW5ub3QgZmluZCBtb2R1bGUgJ1wiK28rXCInXCIpO3Rocm93IGYuY29kZT1cIk1PRFVMRV9OT1RfRk9VTkRcIixmfXZhciBsPW5bb109e2V4cG9ydHM6e319O3Rbb11bMF0uY2FsbChsLmV4cG9ydHMsZnVuY3Rpb24oZSl7dmFyIG49dFtvXVsxXVtlXTtyZXR1cm4gcyhuP246ZSl9LGwsbC5leHBvcnRzLGUsdCxuLHIpfXJldHVybiBuW29dLmV4cG9ydHN9dmFyIGk9dHlwZW9mIHJlcXVpcmU9PVwiZnVuY3Rpb25cIiYmcmVxdWlyZTtmb3IodmFyIG89MDtvPHIubGVuZ3RoO28rKylzKHJbb10pO3JldHVybiBzfSkiLCJtb2R1bGUuZXhwb3J0cyA9IGZ1bmN0aW9uKCkge1xuICAgICAgICAgICAgJF9fcGxhY2Vob2xkZXJfXzBcbiAgICAgICAgICB9LmNhbGwoJF9fcGxhY2Vob2xkZXJfXzEpOyIsIm1vZHVsZS5leHBvcnRzID0ge1xuICAvLyBidWZmZXI6IHJlcXVpcmUoJy4vYnVmZmVyJyksXG4gIC8vIHJtczogcmVxdWlyZSgnLi9ybXMnKSxcbiAgLy8gZW5lcmd5OiByZXF1aXJlKCcuL2VuZXJneScpLFxuICAvLyBjb21wbGV4U3BlY3RydW06IHJlcXVpcmUoJy4vY29tcGxleFNwZWN0cnVtJyksXG4gIC8vIHNwZWN0cmFsU2xvcGU6IHJlcXVpcmUoJy4vc3BlY3RyYWxTbG9wZScpLFxuICAvLyBzcGVjdHJhbENlbnRyb2lkOiByZXF1aXJlKCcuL3NwZWN0cmFsQ2VudHJvaWQnKSxcbiAgLy8gc3BlY3RyYWxSb2xsb2ZmOiByZXF1aXJlKCcuL3NwZWN0cmFsUm9sbG9mZicpLFxuICAvLyBzcGVjdHJhbEZsYXRuZXNzOiByZXF1aXJlKCcuL3NwZWN0cmFsRmxhdG5lc3MnKSxcbiAgLy8gc3BlY3RyYWxTcHJlYWQ6IHJlcXVpcmUoJy4vc3BlY3RyYWxTcHJlYWQnKSxcbiAgLy8gc3BlY3RyYWxTa2V3bmVzczogcmVxdWlyZSgnLi9zcGVjdHJhbFNrZXduZXNzJyksXG4gIC8vIHNwZWN0cmFsS3VydG9zaXM6IHJlcXVpcmUoJy4vc3BlY3RyYWxLdXJ0b3NpcycpLFxuICAvLyBhbXBsaXR1ZGVTcGVjdHJ1bTogcmVxdWlyZSgnLi9hbXBsaXR1ZGVTcGVjdHJ1bScpLFxuICAvLyB6Y3I6IHJlcXVpcmUoJy4vemNyJyksXG4gIC8vIHBvd2VyU3BlY3RydW06IHJlcXVpcmUoJy4vcG93ZXJTcGVjdHJ1bScpLFxuICBsb3VkbmVzczogcmVxdWlyZSgnLi9sb3VkbmVzcycpLFxuICAvLyBwZXJjZXB0dWFsU3ByZWFkOiByZXF1aXJlKCcuL3BlcmNlcHR1YWxTcHJlYWQnKSxcbiAgLy8gcGVyY2VwdHVhbFNoYXJwbmVzczogcmVxdWlyZSgnLi9wZXJjZXB0dWFsU2hhcnBuZXNzJyksXG4gIC8vIG1mY2M6IHJlcXVpcmUoJy4vbWZjYycpXG59OyIsIlxuZnVuY3Rpb24gTG91ZG5lc3Mob3B0cykge1xuICBpZiAoISh0aGlzIGluc3RhbmNlb2YgTG91ZG5lc3MpKSByZXR1cm4gbmV3IExvdWRuZXNzKG9wdHMpO1xuXG4gIC8vIGZlYXR1cmVJbmZvIGZvciBsb3VkbmVzc1xuICB0aGlzLmluZm8gPSB7XG4gICAgXCJ0eXBlXCI6IFwibXVsdGlwbGVBcnJheXNcIixcbiAgICBcImFycmF5TmFtZXNcIjoge1xuICAgICAgXCIxXCI6IFwidG90YWxcIixcbiAgICAgIFwiMlwiOiBcInNwZWNpZmljXCJcbiAgICB9XG4gIH07XG5cbiAgdGhpcy5OVU1fQkFSS19CQU5EUyA9IG9wdHMuTlVNX0JBUktfQkFORFMgfHwgMjQ7XG4gIHRoaXMubm9ybWFsaXNlZFNwZWN0cnVtID0gb3B0cy5ub3JtYWxpc2VkU3BlY3RydW07XG4gIHRoaXMuc2FtcGxlUmF0ZSA9IG9wdHMuc2FtcGxlUmF0ZTtcblxuICB0aGlzLnNwZWNpZmljID0gbmV3IEZsb2F0MzJBcnJheSh0aGlzLk5VTV9CQVJLX0JBTkRTKTtcbiAgdGhpcy5wcm9jZXNzID0gdGhpcy5wcm9jZXNzLmJpbmQodGhpcyk7IC8vIG9wdGltaXplIGxhdGVyXG4gIHRoaXMuYmFya1NjYWxlID0gb3B0cy5iYXJrU2NhbGU7XG4gIHRoaXMuYmJMaW1pdHMgID0gdGhpcy5jb21wdXRlQmFya0JhbmRMaW1pdHModGhpcy5iYXJrU2NhbGUsIHRoaXMubm9ybWFsaXNlZFNwZWN0cnVtLmxlbmd0aCwgdGhpcy5OVU1fQkFSS19CQU5EUyk7XG59XG5cbkxvdWRuZXNzLnByb3RvdHlwZS5jb21wdXRlQmFya0JhbmRMaW1pdHMgPSBmdW5jdGlvbihiYXJrU2NhbGUsIG5TcGVjdHJ1bUxlbmd0aCwgbnVtX2JhcmtfYmFuZHMpIHtcbiAgYmFya1NjYWxlICAgICAgID0gYmFya1NjYWxlICAgICAgIHx8IHRoaXMuYmFya1NjYWxlO1xuICBuU3BlY3RydW1MZW5ndGggPSBuU3BlY3RydW1MZW5ndGggfHwgdGhpcy5ub3JtYWxpc2VkU3BlY3RydW0ubGVuZ3RoO1xuICBudW1fYmFya19iYW5kcyAgPSBudW1fYmFya19iYW5kcyAgfHwgdGhpcy5OVU1fQkFSS19CQU5EUztcbiAgXG4gIHZhciBjdXJyZW50QmFuZEVuZCA9IGJhcmtTY2FsZVtuU3BlY3RydW1MZW5ndGgtMV0vbnVtX2JhcmtfYmFuZHM7XG4gIHZhciBjdXJyZW50QmFuZCA9IDE7XG5cbiAgdmFyIGJiTGltaXRzID0gbmV3IEludDMyQXJyYXkobnVtX2JhcmtfYmFuZHMrMSk7XG4gIGJiTGltaXRzWzBdID0gMDtcblxuICBmb3IodmFyIGkgPSAwOyBpPG5TcGVjdHJ1bUxlbmd0aDsgaSsrKXtcbiAgIHdoaWxlKGJhcmtTY2FsZVtpXSA+IGN1cnJlbnRCYW5kRW5kKSB7XG4gICAgIGJiTGltaXRzW2N1cnJlbnRCYW5kKytdID0gaTtcbiAgICAgY3VycmVudEJhbmRFbmQgPSBjdXJyZW50QmFuZCpiYXJrU2NhbGVbblNwZWN0cnVtTGVuZ3RoLTFdL251bV9iYXJrX2JhbmRzO1xuICAgfVxuICB9XG5cbiAgYmJMaW1pdHNbbnVtX2JhcmtfYmFuZHNdID0gblNwZWN0cnVtTGVuZ3RoLTE7XG5cbiAgcmV0dXJuIGJiTGltaXRzO1xufTtcblxuTG91ZG5lc3MucHJvdG90eXBlLmNvbXB1dGVTcGVjaWZpY0xvdWRuZXNzID0gZnVuY3Rpb24oYmJMaW1pdHMsIHNwZWNpZmljLCBuU3BlY3RydW0sIG51bV9iYXJrX2JhbmRzKSB7XG4gIGJiTGltaXRzICAgICAgICA9IGJiTGltaXRzICAgICAgIHx8wqB0aGlzLmJiTGltaXRzO1xuICBzcGVjaWZpYyAgICAgICAgPSBzcGVjaWZpYyAgICAgICB8fMKgdGhpcy5zcGVjaWZpYztcbiAgblNwZWN0cnVtICAgICAgID0gblNwZWN0cnVtICAgICAgfHwgdGhpcy5ub3JtYWxpc2VkU3BlY3RydW07XG4gIG51bV9iYXJrX2JhbmRzICA9IG51bV9iYXJrX2JhbmRzIHx8IHRoaXMuTlVNX0JBUktfQkFORFM7XG5cbiAgLy8gY29uc29sZS5sb2coYmJMaW1pdHMsIHNwZWNpZmljLCBuU3BlY3RydW0sIG51bV9iYXJrX2JhbmRzKTtcblxuICBmb3IgKHZhciBpID0gMDsgaSA8IG51bV9iYXJrX2JhbmRzOyBpKyspe1xuICAgIHZhciBzdW0gPSAwO1xuXG4gICAgZm9yICh2YXIgaiA9IGJiTGltaXRzW2ldIDsgaiA8IGJiTGltaXRzW2krMV0gOyBqKyspIHtcbiAgICAgIHN1bSArPSBuU3BlY3RydW1bal07XG4gICAgfVxuXG4gICAgc3BlY2lmaWNbaV0gPSBNYXRoLnBvdyhzdW0sIDAuMjMpO1xuICB9XG5cbiAgcmV0dXJuIHNwZWNpZmljO1xufTtcblxuLy8gYmVsb25ncyBpbiB1dGlsc1xuTG91ZG5lc3MucHJvdG90eXBlLnN1bUFycmF5ID0gZnVuY3Rpb24oYXJyYXkpIHtcbiAgdmFyIHN1bSA9IDA7XG5cbiAgZm9yICh2YXIgaSA9IDA7IGkgPCBhcnJheS5sZW5ndGg7IGkrKykge1xuICAgIHN1bSArPSBhcnJheVtpXTtcbiAgfVxuXG4gIHJldHVybiBzdW07XG59O1xuXG5Mb3VkbmVzcy5wcm90b3R5cGUucHJvY2VzcyA9IGZ1bmN0aW9uKCkge1xuICB2YXIgYmFya1NjYWxlID0gdGhpcy5iYXJrU2NhbGU7XG4gIHZhciBuU3BlY3RydW0gPSB0aGlzLm5vcm1hbGlzZWRTcGVjdHJ1bTtcbiAgdmFyIG5TcGVjdHJ1bUxlbmd0aCA9IG5TcGVjdHJ1bS5sZW5ndGg7XG4gIHZhciBOVU1fQkFSS19CQU5EUyA9IHRoaXMuTlVNX0JBUktfQkFORFM7XG4gIHZhciBzcGVjaWZpYyA9IHRoaXMuc3BlY2lmaWM7XG4gIHZhciBiYkxpbWl0cyAgPSB0aGlzLmJiTGltaXRzO1xuXG4gIC8vcHJvY2Vzc1xuICB2YXIgc3BlYyAgPSB0aGlzLmNvbXB1dGVTcGVjaWZpY0xvdWRuZXNzKGJiTGltaXRzLCBzcGVjaWZpYywgblNwZWN0cnVtLCBOVU1fQkFSS19CQU5EUyk7XG4gIC8vIGNvbnNvbGUubG9nKHNwZWMsIHRoaXMuc3BlY2lmaWMpO1xuICB2YXIgdG90YWwgPSB0aGlzLnN1bUFycmF5KHNwZWMpO1xuXG4gIHJldHVybiB7XG4gICAgc3BlY2lmaWM6IHNwZWMsXG4gICAgdG90YWw6IHRvdGFsXG4gIH07XG59O1xuXG5tb2R1bGUuZXhwb3J0cyA9IExvdWRuZXNzOyIsIm1vZHVsZS5leHBvcnR0cyA9IHtcblxuICBcImJ1ZmZlclwiOiB7XG4gICAgXCJ0eXBlXCI6IFwiYXJyYXlcIlxuICB9LFxuICBcInJtc1wiOiB7XG4gICAgXCJ0eXBlXCI6IFwibnVtYmVyXCJcbiAgfSxcbiAgXCJlbmVyZ3lcIjoge1xuICAgIFwidHlwZVwiOiBcIm51bWJlclwiXG4gIH0sXG4gIFwiemNyXCI6IHtcbiAgICBcInR5cGVcIjogXCJudW1iZXJcIlxuICB9LFxuICBcImNvbXBsZXhTcGVjdHJ1bVwiOiB7XG4gICAgXCJ0eXBlXCI6IFwibXVsdGlwbGVBcnJheXNcIixcbiAgICBcImFycmF5TmFtZXNcIjoge1xuICAgICAgXCIxXCI6IFwicmVhbFwiLFxuICAgICAgXCIyXCI6IFwiaW1hZ1wiXG4gICAgfVxuICB9LFxuICBcImFtcGxpdHVkZVNwZWN0cnVtXCI6IHtcbiAgICBcInR5cGVcIjogXCJhcnJheVwiXG4gIH0sXG4gIFwicG93ZXJTcGVjdHJ1bVwiOiB7XG4gICAgXCJ0eXBlXCI6IFwiYXJyYXlcIlxuICB9LFxuICBcInNwZWN0cmFsQ2VudHJvaWRcIjoge1xuICAgIFwidHlwZVwiOiBcIm51bWJlclwiXG4gIH0sXG4gIFwic3BlY3RyYWxGbGF0bmVzc1wiOiB7XG4gICAgXCJ0eXBlXCI6IFwibnVtYmVyXCJcbiAgfSxcbiAgXCJzcGVjdHJhbFNsb3BlXCI6IHtcbiAgICBcInR5cGVcIjogXCJudW1iZXJcIlxuICB9LFxuICBcInNwZWN0cmFsUm9sbG9mZlwiOiB7XG4gICAgXCJ0eXBlXCI6IFwibnVtYmVyXCJcbiAgfSxcbiAgXCJzcGVjdHJhbFNwcmVhZFwiOiB7XG4gICAgXCJ0eXBlXCI6IFwibnVtYmVyXCJcbiAgfSxcbiAgXCJzcGVjdHJhbFNrZXduZXNzXCI6IHtcbiAgICBcInR5cGVcIjogXCJudW1iZXJcIlxuICB9LFxuICBcInNwZWN0cmFsS3VydG9zaXNcIjoge1xuICAgIFwidHlwZVwiOiBcIm51bWJlclwiXG4gIH0sXG4gIFwibG91ZG5lc3NcIjoge1xuICAgIFwidHlwZVwiOiBcIm11bHRpcGxlQXJyYXlzXCIsXG4gICAgXCJhcnJheU5hbWVzXCI6IHtcbiAgICAgIFwiMVwiOiBcInRvdGFsXCIsXG4gICAgICBcIjJcIjogXCJzcGVjaWZpY1wiXG4gICAgfVxuICB9LFxuICBcInBlcmNlcHR1YWxTcHJlYWRcIjoge1xuICAgIFwidHlwZVwiOiBcIm51bWJlclwiXG4gIH0sXG4gIFwicGVyY2VwdHVhbFNoYXJwbmVzc1wiOiB7XG4gICAgXCJ0eXBlXCI6IFwibnVtYmVyXCJcbiAgfSxcbiAgXCJtZmNjXCI6IHtcbiAgICBcInR5cGVcIjogXCJhcnJheVwiXG4gIH1cbn07IiwiLy8gTWV5ZGEgSmF2YXNjcmlwdCBEU1AgbGlicmFyeVxuXG4vLyBEZXBlbmRlbmNpZXNcbi8vIC0tLS0tLS0tLS0tLVxuXG52YXIge0NvbXBsZXhBcnJheX0gXHRcdD0gcmVxdWlyZSgnLi4vbGliL2pzZmZ0L2NvbXBsZXhfYXJyYXknKSxcblx0XHRmZnQgXHRcdFx0XHRcdFx0XHQ9IHJlcXVpcmUoJy4uL2xpYi9qc2ZmdC9mZnQnKTsgLy8gZGVjb3JhdGVzIENvbXBsZXhBcnJheVxuXG52YXIge2lzUG93ZXJPZlR3bywgwrV9ID0gcmVxdWlyZSgnLi91dGlscycpO1xuXG5cbi8vIENvbnN0cnVjdG9yXG4vLyAtLS0tLS0tLS0tLVxuXG5jbGFzcyBNZXlkYSB7XG5cblx0Y29uc3RydWN0b3IoYXVkaW9Db250ZXh0LCBzcmMsIGJ1ZlNpemUsIGNhbGxiYWNrKSB7XG5cdFx0XG5cdFx0Ly8gbWFrZSBzdXJlIHdlIGNhbiB3b3JrXG5cdFx0aWYoIWlzUG93ZXJPZlR3byhidWZTaXplKSkge1xuXHRcdFx0dGhyb3cgbmV3IEVycm9yKFwiQnVmZmVyIHNpemUgaXMgbm90IGEgcG93ZXIgb2YgdHdvOiBNZXlkYSB3aWxsIG5vdCBydW4uXCIpO1xuXHRcdH1cblx0XHRcblx0XHRpZighYXVkaW9Db250ZXh0KSB7XG5cdFx0XHR0aHJvdyBuZXcgRXJyb3IoXCJBdWRpb0NvbnRleHQgd2Fzbid0IHNwZWNpZmllZDogTWV5ZGEgd2lsbCBub3QgcnVuLlwiKTtcblx0XHR9XG5cblx0XHR2YXIgYnVmZmVyU2l6ZSA9IGJ1ZlNpemUgfHwgMjU2OyAvL2RlZmF1bHQgYnVmZmVyIHNpemVcblx0XHR2YXIgc2FtcGxlUmF0ZSA9IGF1ZGlvQ29udGV4dC5zYW1wbGVSYXRlO1xuXHRcdHZhciBzb3VyY2UgPSBzcmM7IC8vaW5pdGlhbCBzb3VyY2VcblxuXHRcdHRoaXMuYXVkaW9Db250ZXh0ID0gYXVkaW9Db250ZXh0O1xuXHRcdHRoaXMuZmVhdHVyZUV4dHJhY3RvcnMgPSB7fTtcblxuXHRcdC8vY2FsbGJhY2sgY29udHJvbGxlcnNcblx0XHR0aGlzLkVYVFJBQ1RJT05fU1RBUlRFRCA9IGZhbHNlO1xuXHRcdHRoaXMuX2ZlYXR1cmVzVG9FeHRyYWN0ID0gbnVsbDtcblx0XHRcblx0XHQvL1dJTkRPV0lOR1xuXHRcdC8vc2V0IGRlZmF1bHRcblx0XHR0aGlzLndpbmRvd2luZ0Z1bmN0aW9uID0gXCJoYW5uaW5nXCI7XG5cblx0XHQvL2luaXRpbGl6ZSBiYXJrIHNjYWxlICh0byBiZSB1c2VkIGluIG1vc3QgcGVyY2VwdHVhbCBmZWF0dXJlcykuXG5cdFx0dGhpcy5iYXJrU2NhbGUgPSB0aGlzLmNvbXB1dGVCYXJrU2NhbGUoYnVmZmVyU2l6ZSwgc2FtcGxlUmF0ZSk7XG5cblx0XHQvL2NyZWF0ZSB3aW5kb3dzXG5cdFx0dGhpcy5oYW5uaW5nID0gdGhpcy5jb21wdXRlSGFubmluZyhidWZmZXJTaXplKTtcblx0XHR0aGlzLmhhbW1pbmcgPSB0aGlzLmNvbXB1dGVIYW1taW5nKGJ1ZmZlclNpemUpO1xuXHRcdC8vIHRoaXMuYmxhY2ttYW4gPSB0aGlzLmNvbXB1dGVCbGFja21hbihidWZmZXJTaXplKTtcblxuXHRcdHRoaXMuZmVhdHVyZUluZm8gPSByZXF1aXJlKCcuL2ZlYXR1cmUtaW5mbycpOyAvLyB0byBiZSBpbmNsdWRlZCBwZXIgbW9kdWxlXG5cblx0XHQvL2NyZWF0ZSBjb21wbGV4YXJyYXkgdG8gaG9sZCB0aGUgc3BlY3RydW1cblx0XHR2YXIgY29tcHV0ZWRTcGVjdHJ1bURhdGEgPSB0aGlzLmNvbXB1dGVTcGVjdHJ1bURhdGEoYnVmZmVyU2l6ZSk7XG5cblx0XHQvL2Fzc2lnbiB0byBtZXlkYVxuXHRcdHRoaXMuc3BlY3RydW1EYXRhID0gY29tcHV0ZWRTcGVjdHJ1bURhdGEuZGF0YTtcblx0XHR0aGlzLmNvbXBsZXhTcGVjdHJ1bSA9IGNvbXB1dGVkU3BlY3RydW1EYXRhLnNwZWN0cnVtO1xuXHRcdHRoaXMuYW1wU3BlY3RydW0gPSBjb21wdXRlZFNwZWN0cnVtRGF0YS5hbXBTcGVjdHJ1bTtcblxuXHRcdC8vIGNvbnNvbGUubG9nKGNvbXB1dGVkU3BlY3RydW1EYXRhKTtcblxuXHRcdFx0Ly8gY29uc29sZS5sb2coYW1wU3BlY3RydW0pXG5cdFx0dGhpcy5pbml0aWFsaXNlRXh0cmFjdG9ycygpO1xuXG5cdFx0Ly9jcmVhdGUgbm9kZXNcblx0XHR3aW5kb3cuc3BuID0gYXVkaW9Db250ZXh0LmNyZWF0ZVNjcmlwdFByb2Nlc3NvcihidWZmZXJTaXplLDEsMSk7XG5cblx0XHR3aW5kb3cuc3BuLm9uYXVkaW9wcm9jZXNzID0gKGUpID0+IHtcblxuXHRcdFx0Ly90aGlzIGlzIHRvIG9idGFpbiB0aGUgY3VycmVudCBhbXBsaXR1ZGUgc3BlY3RydW1cblx0XHRcdHZhciBzaWduYWwgPSB0aGlzLnNpZ25hbCA9IGUuaW5wdXRCdWZmZXIuZ2V0Q2hhbm5lbERhdGEoMCk7XG5cdFx0XHR2YXIgZGF0YSA9IHRoaXMuc3BlY3RydW1EYXRhO1xuXHRcdFx0dmFyIHNwZWMgPSB0aGlzLmNvbXBsZXhTcGVjdHJ1bTtcblx0XHRcdHZhciBhbXBTcGVjdHJ1bSA9IHRoaXMuYW1wU3BlY3RydW07XG5cdFx0XHR2YXIgd2luZG93ZWRTaWduYWwgPSB0aGlzLmNvbXB1dGVXaW5kb3coc2lnbmFsLCB0aGlzLndpbmRvd2luZ0Z1bmN0aW9uKTtcblxuXHRcdFx0Ly9tYXAgdGltZSBkb21haW5cblx0XHRcdGRhdGEubWFwKGZ1bmN0aW9uKHZhbHVlLCBpLCBuKSB7XG5cdFx0XHRcdHZhbHVlLnJlYWwgPSB3aW5kb3dlZFNpZ25hbFtpXTtcblx0XHRcdH0pO1xuXG5cdFx0XHQvL2NhbGN1bGF0ZSBhbXBsaXR1ZGVcblx0XHRcdHRoaXMuY29tcHV0ZUFtcGxpdHVkZShzcGVjLCBhbXBTcGVjdHJ1bSwgYnVmZmVyU2l6ZSk7XG5cblx0XHRcdC8vY2FsbCBjYWxsYmFjayBpZiBhcHBsaWNhYmxlXG5cdFx0XHRpZiAodHlwZW9mIGNhbGxiYWNrID09PSBcImZ1bmN0aW9uXCIgJiYgRVhUUkFDVElPTl9TVEFSVEVEKSB7XG5cdFx0XHRcdGNhbGxiYWNrKHRoaXMuZ2V0KHRoaXMuX2ZlYXR1cmVzVG9FeHRyYWN0KSk7XG5cdFx0XHR9XG5cblx0XHR9O1xuXG5cdFx0d2luZG93LnNwbi5jb25uZWN0KGF1ZGlvQ29udGV4dC5kZXN0aW5hdGlvbik7XG5cdFx0c291cmNlLmNvbm5lY3Qod2luZG93LnNwbiwgMCwgMCk7XG5cdFx0XG5cdFx0Ly8gY29uc3RydWN0b3JzIHJldHVybiBcInRoaXNcIiBieSBkZWZhdWx0XG5cdH1cblxuXG5cblx0Ly8gQ29tcHV0ZSBtZXRob2RzXG5cdC8vIC0tLS0tLS0tLS0tLS0tLVxuXG5cdGNvbXB1dGVBbXBsaXR1ZGUoY29tcGxleFNwZWN0cnVtLCBhbXBTcGVjdHJ1bSwgYnVmZmVyU2l6ZSkge1xuXG5cdFx0Ly8gd29ya3MgaW5zaWRlIGFuZCBvdXRzaWRlIE1leWRhXG5cdFx0Y29tcGxleFNwZWN0cnVtID0gY29tcGxleFNwZWN0cnVtIHx8wqB0aGlzLmNvbXBsZXhTcGVjdHJ1bTtcblx0XHRhbXBTcGVjdHJ1bSBcdFx0PSBhbXBTcGVjdHJ1bSBcdFx0fHzCoHRoaXMuYW1wU3BlY3RydW07XG5cdFx0YnVmZmVyU2l6ZSBcdFx0XHQ9IGJ1ZmZlclNpemUgXHRcdFx0fHzCoHRoaXMuYnVmZmVyU2l6ZTtcblxuXHRcdGZvciAodmFyIGkgPSAwOyBpIDwgYnVmZmVyU2l6ZS8yOyBpKyspIHtcblx0XHRcdGFtcFNwZWN0cnVtW2ldID0gTWF0aC5zcXJ0KE1hdGgucG93KGNvbXBsZXhTcGVjdHJ1bS5yZWFsW2ldLDIpICsgTWF0aC5wb3coY29tcGxleFNwZWN0cnVtLmltYWdbaV0sMikpO1xuXHRcdH1cblx0fVxuXG5cdGNvbXB1dGVIYW1taW5nKGJ1ZmZlclNpemUpIHtcblx0XHRidWZmZXJTaXplID0gYnVmZmVyU2l6ZSB8fMKgdGhpcy5idWZmZXJTaXplO1xuXG5cdFx0dmFyIGhhbW1pbmcgPSBuZXcgRmxvYXQzMkFycmF5KGJ1ZmZlclNpemUpO1xuXHRcdGZvciAodmFyIGkgPSAwOyBpIDwgYnVmZmVyU2l6ZTsgaSsrKSB7XG5cdFx0XHQvL0FjY29yZGluZyB0byBodHRwOi8vdWsubWF0aHdvcmtzLmNvbS9oZWxwL3NpZ25hbC9yZWYvaGFtbWluZy5odG1sXG5cdFx0XHRoYW1taW5nW2ldID0gMC41NCAtIDAuNDYqTWF0aC5jb3MoMipNYXRoLlBJKihpL2J1ZmZlclNpemUtMSkpO1xuXHRcdH1cblxuXHRcdHJldHVybiBoYW1taW5nO1xuXHRcdH1cblxuXHRjb21wdXRlSGFubmluZyhidWZmZXJTaXplKSB7XG5cdFx0YnVmZmVyU2l6ZSA9IGJ1ZmZlclNpemUgfHzCoHRoaXMuYnVmZmVyU2l6ZTtcblxuXHRcdHZhciBoYW5uaW5nID0gbmV3IEZsb2F0MzJBcnJheShidWZmZXJTaXplKTtcblx0XHRmb3IgKHZhciBpID0gMDsgaSA8IGJ1ZmZlclNpemU7IGkrKykge1xuXHRcdFx0Ly9BY2NvcmRpbmcgdG8gdGhlIFIgZG9jdW1lbnRhdGlvbiBodHRwOi8vcmdtLm9nYWxhYi5uZXQvUkdNL1JfcmRmaWxlP2Y9R0VORUFyZWFkL21hbi9oYW5uaW5nLndpbmRvdy5SZCZkPVJfQ0Ncblx0XHRcdGhhbm5pbmdbaV0gPSAwLjUgLSAwLjUqTWF0aC5jb3MoMipNYXRoLlBJKmkvKGJ1ZmZlclNpemUtMSkpO1xuXHRcdH1cblxuXHRcdHJldHVybiBoYW5uaW5nO1xuXHR9XG5cblx0Ly9VTkZJTklTSEVEIC0gYmxhY2ttYW4gd2luZG93IGltcGxlbWVudGF0aW9uXG5cdC8qXG5cdFx0Y29tcHV0ZUJsYWNrbWFuKGJ1ZmZlclNpemUpIHtcblx0XHRidWZmZXJTaXplID0gYnVmZmVyU2l6ZSB8fMKgdGhpcy5idWZmZXJTaXplO1xuXHRcdFxuXHRcdHZhciBibGFja21hbiA9IG5ldyBGbG9hdDMyQXJyYXkoYnVmZmVyU2l6ZSk7XG5cdFx0Ly9BY2NvcmRpbmcgdG8gaHR0cDovL3VrLm1hdGh3b3Jrcy5jb20vaGVscC9zaWduYWwvcmVmL2JsYWNrbWFuLmh0bWxcblx0XHQvL2ZpcnN0IGhhbGYgb2YgdGhlIHdpbmRvd1xuXHRcdGZvciAodmFyIGkgPSAwOyBpIDwgKGJ1ZmZlclNpemUgJSAyKSA/IChidWZmZXJTaXplKzEpLzIgOiBidWZmZXJTaXplLzI7IGkrKykge1xuXHRcdFx0dGhpcy5ibGFja21hbltpXSA9IDAuNDIgLSAwLjUqTWF0aC5jb3MoMipNYXRoLlBJKmkvKGJ1ZmZlclNpemUtMSkpICsgMC4wOCpNYXRoLmNvcyg0Kk1hdGguUEkqaS8oYnVmZmVyU2l6ZS0xKSk7XG5cdFx0fVxuXHRcdC8vc2Vjb25kIGhhbGYgb2YgdGhlIHdpbmRvd1xuXHRcdGZvciAodmFyIGkgPSBidWZmZXJTaXplLzI7IGkgPiAwOyBpLS0pIHtcblx0XHRcdHRoaXMuYmxhY2ttYW5bYnVmZmVyU2l6ZSAtIGldID0gdGhpcy5ibGFja21hbltpXTtcblx0XHR9XG5cdH07XG5cdCovXG5cblx0Y29tcHV0ZVdpbmRvdyhzaWcsIHR5cGUpIHtcblxuXHRcdHZhciBpLCBsZW4gPSBzaWcubGVuZ3RoO1xuXHRcdHZhciB3aW5kb3dlZCA9IG5ldyBGbG9hdDMyQXJyYXkobGVuKTtcblxuXHRcdGZvciAoaSA9IDA7IGkgPCBsZW47IGkrKykge1xuXHRcdFx0d2luZG93ZWRbaV0gPSBzaWdbaV0gKiB0aGlzW3R5cGVdW2ldO1xuXHRcdH1cblx0XHRcblx0XHRyZXR1cm4gd2luZG93ZWQ7XG5cdFx0fVxuXG5cdFx0Y29tcHV0ZUJhcmtTY2FsZShidWZmZXJTaXplLCBzYW1wbGVSYXRlKSB7XG5cdFx0YnVmZmVyU2l6ZSA9IGJ1ZmZlclNpemUgfHzCoHRoaXMuYnVmZmVyU2l6ZTtcblx0XHRzYW1wbGVSYXRlID0gc2FtcGxlUmF0ZSB8fMKgdGhpcy5zYW1wbGVSYXRlO1xuXG5cdCAgdmFyIGJhcmtTY2FsZSA9IG5ldyBGbG9hdDMyQXJyYXkoYnVmZmVyU2l6ZSk7XG5cblx0ICBmb3IodmFyIGkgPSAwOyBpIDwgYnVmZmVyU2l6ZTsgaSsrKXtcblx0ICAgIGJhcmtTY2FsZVtpXSA9IGkgKiBzYW1wbGVSYXRlIC8gKGJ1ZmZlclNpemUpO1xuXHQgICAgYmFya1NjYWxlW2ldID0gMTMgKiBNYXRoLmF0YW4oYmFya1NjYWxlW2ldLzEzMTUuOCkgKyAzLjUgKiBNYXRoLmF0YW4oTWF0aC5wb3coKGJhcmtTY2FsZVtpXS83NTE4KSwgMikpO1xuXHQgIH1cblxuXHQgIHJldHVybiBiYXJrU2NhbGU7XG5cdH1cblxuXHRjb21wdXRlU3BlY3RydW1EYXRhKGJ1ZmZlclNpemUpIHtcblx0XHRidWZmZXJTaXplID0gYnVmZmVyU2l6ZSB8fMKgdGhpcy5idWZmZXJTaXplO1xuXG5cdFx0Ly9jcmVhdGUgY29tcGxleGFycmF5IHRvIGhvbGQgdGhlIHNwZWN0cnVtXG5cdFx0dmFyIGRhdGEgPSBuZXcgQ29tcGxleEFycmF5KGJ1ZmZlclNpemUpO1xuXHRcdHZhciBzcGVjdHJ1bSA9IGRhdGEuRkZUKCk7IC8vdHJhbnNmb3JtXG5cdFx0dmFyIGFtcFNwZWN0cnVtID0gbmV3IEZsb2F0MzJBcnJheShidWZmZXJTaXplLzIpO1xuXG5cdFx0cmV0dXJuIHtcblx0XHRcdGRhdGE6IGRhdGEsXG5cdFx0XHRzcGVjdHJ1bTogc3BlY3RydW0sXG5cdFx0XHRhbXBTcGVjdHJ1bTogYW1wU3BlY3RydW1cblx0XHR9O1xuXHR9XG5cblxuXHQvLyBNZXlkYSBtZXRob2RzXG5cdC8vIC0tLS0tLS0tLS0tLS1cblx0XG5cdC8vIGxvYWRzIGFsbCB0aGUgZXh0cmFjdG9yIG9iamVjdHNcblx0Ly8gaW5pdGlhbGl6ZXMgdGhlbSBhbmQgYmluZHMgdGhlbVxuXHQvLyB0byB0aGUgZmVhdHVyZUV4dHJhY3RvcnNcblx0Ly8gYW5kIGZlYXR1cmVJbmZvIGxpc3RzXG5cdC8vIEBOT1RFIChjL3NoKW91bGQgYmUgaGFuZGVsZWQgZGlmZmVyZW50bHlcblx0aW5pdGlhbGlzZUV4dHJhY3RvcnMoKSB7XG5cblx0XHR2YXIgZXh0cmFjdG9ycyA9IHJlcXVpcmUoJy4vZXh0cmFjdG9ycycpO1xuXG5cdFx0Ly8gTG91ZG5lc3Ncblx0XHR2YXIgbG91ZG5lc3MgPSBleHRyYWN0b3JzLmxvdWRuZXNzKHtcblx0XHRcdE5VTV9CQVJLX0JBTkRTOiAyNCxcblx0XHRcdGJhcmtTY2FsZTogdGhpcy5iYXJrU2NhbGUsXG5cdFx0XHRub3JtYWxpc2VkU3BlY3RydW06IHRoaXMuYW1wU3BlY3RydW0sXG5cdFx0XHRzYW1wbGVSYXRlOiB0aGlzLmF1ZGlvQ29udGV4dC5zYW1wbGVSYXRlXG5cdFx0fSk7XG5cblx0XHR0aGlzLmZlYXR1cmVFeHRyYWN0b3JzLmxvdWRuZXNzID0gbG91ZG5lc3M7XG5cdFx0dGhpcy5mZWF0dXJlSW5mby5sb3VkbmVzcyA9IGxvdWRuZXNzLmluZm87XG5cblx0XHQvLyBSZXN0IG9mIHRoZSBleHRyYWN0b3JzXG5cdFx0Ly8gLi4uXG5cdH1cblxuXG5cdC8vc291cmNlIHNldHRlciBtZXRob2Rcblx0c2V0U291cmNlKF9zcmMpIHtcblx0XHRfc3JjLmNvbm5lY3Qod2luZG93LnNwbik7XG5cdH1cblxuXHRzdGFydChmZWF0dXJlcykge1xuXHRcdHRoaXMuX2ZlYXR1cmVzVG9FeHRyYWN0ID0gZmVhdHVyZXM7XG5cdFx0dGhpcy5FWFRSQUNUSU9OX1NUQVJURUQgPSB0cnVlO1xuXHR9XG5cblx0c3RvcCgpIHtcblx0XHR0aGlzLl9mZWF0dXJlc1RvRXh0cmFjdCA9IG51bGw7XG5cdFx0dGhpcy5FWFRSQUNUSU9OX1NUQVJURUQgPSBmYWxzZTtcblx0fVxuXG4vLyBkYXRhIHB1bGxpbmdcblx0Z2V0KGZlYXR1cmUpIHtcblxuXHRcdGlmKHR5cGVvZiBmZWF0dXJlID09PSBcIm9iamVjdFwiKXtcblx0XHRcdHZhciByZXN1bHRzID0ge307XG5cdFx0XHRmb3IgKHZhciB4ID0gMDsgeCA8IGZlYXR1cmUubGVuZ3RoOyB4Kyspe1xuXHRcdFx0XHR0cnl7XG5cdFx0XHRcdFx0cmVzdWx0c1tmZWF0dXJlW3hdXSA9ICh0aGlzLmZlYXR1cmVFeHRyYWN0b3JzW2ZlYXR1cmVbeF1dLnByb2Nlc3ModGhpcy5zaWduYWwpKTtcblx0XHRcdFx0fSBjYXRjaCAoZSl7XG5cdFx0XHRcdFx0Y29uc29sZS5lcnJvcihlKTtcblx0XHRcdFx0fVxuXHRcdFx0fVxuXHRcdFx0cmV0dXJuIHJlc3VsdHM7XG5cdFx0fSBlbHNlIGlmICh0eXBlb2YgZmVhdHVyZSA9PT0gXCJzdHJpbmdcIil7XG5cdFx0XHRyZXR1cm4gdGhpcy5mZWF0dXJlRXh0cmFjdG9yc1tmZWF0dXJlXS5wcm9jZXNzKHRoaXMuc2lnbmFsKTtcblx0XHR9IGVsc2V7XG5cdFx0XHR0aHJvdyBuZXcgRXJyb3IoXCJJbnZhbGlkIEZlYXR1cmUgRm9ybWF0XCIpO1xuXHRcdH1cblx0fVxuXG59XG5cbm1vZHVsZS5leHBvcnRzID0gTWV5ZGE7XG4iLCJtb2R1bGUuZXhwb3J0cy7CtSA9IGZ1bmN0aW9uKGksIGFtcGxpdHVkZVNwZWN0KXtcbiAgdmFyIG51bWVyYXRvciA9IDA7XG4gIHZhciBkZW5vbWluYXRvciA9IDA7XG4gIFxuICBmb3IodmFyIGsgPSAwOyBrIDwgYW1wbGl0dWRlU3BlY3QubGVuZ3RoOyBrKyspe1xuICAgIG51bWVyYXRvciArPSBNYXRoLnBvdyhrLGkpKk1hdGguYWJzKGFtcGxpdHVkZVNwZWN0W2tdKTtcbiAgICBkZW5vbWluYXRvciArPSBhbXBsaXR1ZGVTcGVjdFtrXTtcbiAgfVxuXG4gIHJldHVybiBudW1lcmF0b3IvZGVub21pbmF0b3I7XG59O1xuXG5tb2R1bGUuZXhwb3J0cy5pc1Bvd2VyT2ZUd28gPSBmdW5jdGlvbihudW0pIHtcbiAgd2hpbGUgKCgobnVtICUgMikgPT09IDApICYmIG51bSA+IDEpIHtcbiAgICBudW0gLz0gMjtcbiAgfVxuXG4gIHJldHVybiAobnVtID09IDEpO1xufTsiXX0=

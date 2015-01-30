
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
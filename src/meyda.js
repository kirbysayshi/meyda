// Meyda Javascript DSP library

// Dependencies
// ------------

var {ComplexArray} 		= require('../lib/jsfft/complex_array'),
		fft 							= require('../lib/jsfft/fft'); // decorates ComplexArray

var {isPowerOfTwo, µ} = require('./utils');


// Constructor
// -----------

class Meyda {

	constructor(audioContext, src, bufSize, callback) {
		
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

		window.spn.onaudioprocess = (e) => {

			//this is to obtain the current amplitude spectrum
			var signal = this.signal = e.inputBuffer.getChannelData(0);
			var data = this.spectrumData;
			var spec = this.complexSpectrum;
			var ampSpectrum = this.ampSpectrum;
			var windowedSignal = this.computeWindow(signal, this.windowingFunction);

			//map time domain
			data.map(function(value, i, n) {
				value.real = windowedSignal[i];
			});

			//calculate amplitude
			this.computeAmplitude(spec, ampSpectrum, bufferSize);

			//call callback if applicable
			if (typeof callback === "function" && EXTRACTION_STARTED) {
				callback(this.get(this._featuresToExtract));
			}

		};

		window.spn.connect(audioContext.destination);
		source.connect(window.spn, 0, 0);
		
		// constructors return "this" by default
	}



	// Compute methods
	// ---------------

	computeAmplitude(complexSpectrum, ampSpectrum, bufferSize) {

		// works inside and outside Meyda
		complexSpectrum = complexSpectrum || this.complexSpectrum;
		ampSpectrum 		= ampSpectrum 		|| this.ampSpectrum;
		bufferSize 			= bufferSize 			|| this.bufferSize;

		for (var i = 0; i < bufferSize/2; i++) {
			ampSpectrum[i] = Math.sqrt(Math.pow(complexSpectrum.real[i],2) + Math.pow(complexSpectrum.imag[i],2));
		}
	}

	computeHamming(bufferSize) {
		bufferSize = bufferSize || this.bufferSize;

		var hamming = new Float32Array(bufferSize);
		for (var i = 0; i < bufferSize; i++) {
			//According to http://uk.mathworks.com/help/signal/ref/hamming.html
			hamming[i] = 0.54 - 0.46*Math.cos(2*Math.PI*(i/bufferSize-1));
		}

		return hamming;
		}

	computeHanning(bufferSize) {
		bufferSize = bufferSize || this.bufferSize;

		var hanning = new Float32Array(bufferSize);
		for (var i = 0; i < bufferSize; i++) {
			//According to the R documentation http://rgm.ogalab.net/RGM/R_rdfile?f=GENEAread/man/hanning.window.Rd&d=R_CC
			hanning[i] = 0.5 - 0.5*Math.cos(2*Math.PI*i/(bufferSize-1));
		}

		return hanning;
	}

	//UNFINISHED - blackman window implementation
	/*
		computeBlackman(bufferSize) {
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

	computeWindow(sig, type) {

		var i, len = sig.length;
		var windowed = new Float32Array(len);

		for (i = 0; i < len; i++) {
			windowed[i] = sig[i] * this[type][i];
		}
		
		return windowed;
		}

		computeBarkScale(bufferSize, sampleRate) {
		bufferSize = bufferSize || this.bufferSize;
		sampleRate = sampleRate || this.sampleRate;

	  var barkScale = new Float32Array(bufferSize);

	  for(var i = 0; i < bufferSize; i++){
	    barkScale[i] = i * sampleRate / (bufferSize);
	    barkScale[i] = 13 * Math.atan(barkScale[i]/1315.8) + 3.5 * Math.atan(Math.pow((barkScale[i]/7518), 2));
	  }

	  return barkScale;
	}

	computeSpectrumData(bufferSize) {
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
	}


	// Meyda methods
	// -------------
	
	// loads all the extractor objects
	// initializes them and binds them
	// to the featureExtractors
	// and featureInfo lists
	// @NOTE (c/sh)ould be handeled differently
	initialiseExtractors() {

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
	}


	//source setter method
	setSource(_src) {
		_src.connect(window.spn);
	}

	start(features) {
		this._featuresToExtract = features;
		this.EXTRACTION_STARTED = true;
	}

	stop() {
		this._featuresToExtract = null;
		this.EXTRACTION_STARTED = false;
	}

// data pulling
	get(feature) {

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
	}

}

module.exports = Meyda;

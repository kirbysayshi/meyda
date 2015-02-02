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

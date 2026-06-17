function Pitchshift(fftFrameSize, sampleRate, algo) {
	if (arguments.length)
		this.getready(fftFrameSize, sampleRate, algo);
}

Pitchshift.prototype.getready = function(fftFrameSize, sampleRate, algo) {
	this.fftFrameSize_ = fftFrameSize;
	this.sampleRate_ = sampleRate;
	this.hannWindow_ = []
	this.gRover_ = false;
	this.algo = algo || "FFT";
	//This has to go.
	this.MAX_FRAME_LENGTH = 8192

	function newFilledArray(length, val) {
		var intLength = Math.floor(length);
		var array = [];
		for (var i = 0; i < intLength; i++)
			array[i] = val;
		return array;
	}

	this.gInFIFO = newFilledArray(this.MAX_FRAME_LENGTH, 0);
	this.gOutFIFO = newFilledArray(this.MAX_FRAME_LENGTH, 0);
	this.gLastPhase = newFilledArray(this.MAX_FRAME_LENGTH / 2 + 1, 0);
	this.gSumPhase = newFilledArray(this.MAX_FRAME_LENGTH / 2 + 1, 0);
	this.gOutputAccum = newFilledArray(2 * this.MAX_FRAME_LENGTH, 0);
	this.gAnaFreq = newFilledArray(this.MAX_FRAME_LENGTH, 0);
	this.gAnaMagn = newFilledArray(this.MAX_FRAME_LENGTH, 0);
	this.gSynFreq = newFilledArray(this.MAX_FRAME_LENGTH, 0);
	this.gSynMagn = newFilledArray(this.MAX_FRAME_LENGTH, 0);
	//this.gFFTworksp = newFilledArray(2 * this.MAX_FRAME_LENGTH, 0);
	//Not two, 'cos we haven't to fill phases with 0's.
	this.gFFTworksp = newFilledArray(this.fftFrameSize_, 0);
	//Real and imaginary parts of the resynthesized signal
	this.real_ = [];
	this.imag_ = [];
	//Output data.
	this.outdata = [];
	this.hannWindow_ = [];
	//Pre-generating Hann wavetable
	for (k = 0; k < fftFrameSize; k++)
		this.hannWindow_[k] = WindowFunction.Hann(fftFrameSize, k);
	//Init once, use always.
	if (this.algo === "FFT")
		this.fft = new FFT(this.fftFrameSize_, this.sampleRate_);
	else if (this.algo === "RFFT")
		this.fft = new RFFT(this.fftFrameSize_, this.sampleRate_);
	else
		throw new Error("Invalid DFT algorithm selected " + this.algo);
	//Probably we don't need this.
	this.invFFT = new FFT(this.fftFrameSize_, this.sampleRate_);

	console.log("Pitchshift.prototype.getready returns back");
};

Pitchshift.prototype.process = function(pitchShift, numSampsToProcess, osamp, indata) {
	function setArray(array, length, val) {
		var intLength = Math.floor(length);
		for (var i = 0; i < intLength; i++)
			array[i] = val;
	}
	//PitchShift: factor value which is between 0.5 (one octave down) and 2. (one octave up).
	var fftFrameSize2 = this.fftFrameSize_ / 2, stepSize = this.fftFrameSize_ / osamp,
		freqPerBin = this.sampleRate_ / this.fftFrameSize_, expct = 2. * Math.PI * stepSize / this.fftFrameSize_,
		inFifoLatency = this.fftFrameSize_ - stepSize, j, k = 0, magn, phase, tmp, qpd, index, signal;
	if (this.gRover_ === false)
		this.gRover_ = inFifoLatency;
	//Main processing loop
	for (j = 0; j < numSampsToProcess; j++) {
		//As long as we have not yet collected enough data just read in
		this.gInFIFO[this.gRover_] = indata[j];
		this.outdata[j] = this.gOutFIFO[this.gRover_ - inFifoLatency];
		this.gRover_++;

		//We have enough data for processing
		if (this.gRover_ >= this.fftFrameSize_) {
			this.gRover_ = inFifoLatency;
			//Windowing
			for (k = 0; k < this.fftFrameSize_; k++) {
				//Need the signal for the FFT.
				this.gFFTworksp[k] = this.gInFIFO[k] * this.hannWindow_[k];
				//this.gFFTworksp[k][1] = 0.;
			}
			this.fft.forward(this.gFFTworksp);

			//Analysis
			for (k = 0; k <= fftFrameSize2; k++) {
				//These ifs make the pitchshifter code dependent on the DFT implementation; we should decorate DFTs instead.
				if (this.algo === "FFT") {
					//Taking some "private" member out of fft here.
					magn = 2 * Math.sqrt(this.fft.real[k] * this.fft.real[k] + this.fft.imag[k] * this.fft.imag[k]);
					//aka magn = spectrum[k];
					phase = Math.atan2(this.fft.imag[k], this.fft.real[k]);
				} else if (this.algo === "RFFT") {
					//Because having the same interface but a different output schema
					//in the same library is a great fucking idea!
					var imaginary, real;
					real = this.fft.trans[k];
					imaginary = k==0 ? 0 : this.fft.trans[this.fftFrameSize_ - k];
					magn = 2 * Math.sqrt(real * real + imaginary * imaginary);
					phase = Math.atan2(imaginary, real);
				} else {
					//If we used the constructor, we can't be here.
					throw new Error("Invalid DFT algorithm selected " + this.algo);
				}

				//Compute phase difference
				tmp = phase - this.gLastPhase[k];
				this.gLastPhase[k] = phase;
				//Subtract expected phase difference
				tmp -= k * expct;

				//Map delta phase into +/- Pi interval
				/* Floor and ceil should emulate the behaviour of a C float -> long int conversion
				"Truncating conversion means that any fractional part is discarded, so that e.g.
				3.9 is converted to 3". (http://www.cs.tut.fi/~jkorpela/round.html)*/

				qpd = tmp / Math.PI;
				if (qpd >= 0) {
					qpd = Math.floor(qpd);
					/* This probably won't work like in C */
					qpd += qpd & 1;
				} else {
					qpd = Math.ceil(qpd);
					qpd -= qpd & 1;
				}
				tmp -= Math.PI * qpd;
				//Get deviation from bin frequency from the +/- Pi interval
				tmp = osamp * tmp / (2 * Math.PI);
				//Compute the k-th partials' true frequency
				tmp = k * freqPerBin + tmp * freqPerBin;
				//Store magnitude and true frequency in analysis arrays
				this.gAnaMagn[k] = magn;
				this.gAnaFreq[k] = tmp;
			}

			//Pitch shifting
			//memset(gSynMagn, 0, fftFrameSize*sizeof(float));
			//memset(gSynFreq, 0, fftFrameSize*sizeof(float));
			setArray(this.gSynMagn, this.fftFrameSize_, 0);
			setArray(this.gSynFreq, this.fftFrameSize_, 0);
			for (k = 0; k <= fftFrameSize2; k++) {
				//This is an int multiplication in C.
				index = Math.floor(k * pitchShift);
				if (index <= fftFrameSize2) {
					this.gSynMagn[index] += this.gAnaMagn[k];
					this.gSynFreq[index] = this.gAnaFreq[k] * pitchShift;
				}
			}

			//Synthesis
			for (k = 0; k <= fftFrameSize2; k++) {
				//Get magnitude and true frequency from synthesis arrays
				magn = this.gSynMagn[k];
				tmp = this.gSynFreq[k];
				//Subtract bin mid frequency
				tmp -= k * freqPerBin;
				//Get bin deviation from freq deviation
				tmp /= freqPerBin;
				//Take osamp into account
				tmp = 2. * Math.PI * tmp / osamp;
				//Add the overlap phase advance back in
				tmp += k * expct;
				//Accumulate delta phase to get bin phase
				this.gSumPhase[k] += tmp;
				phase = this.gSumPhase[k];
				//Get real and imag part
				this.real_[k] = magn * Math.cos(phase);
				this.imag_[k] = magn * Math.sin(phase);
			}
			//Zero negative frequencies
			for (k = ((fftFrameSize2) + 1); (k < this.fftFrameSize_); k++) {
				//That's ok, otherwise inverse fft has a fit.
				this.real_[k] = 0;
				this.imag_[k] = 0;
			}
			//Do the Inverse transform
			signal = this.invFFT.inverse(this.real_, this.imag_);
			//Do inverse windowing and add to output accumulator
			for (k = 0; k < this.fftFrameSize_; k++)
				this.gOutputAccum[k] += this.hannWindow_[k] * signal[k];
			for (k = 0; k < stepSize; k++)
				this.gOutFIFO[k] = this.gOutputAccum[k];
			//Shift the output accumulator.
			//Rough memmove implementation.
			var tempArray = this.gOutputAccum.slice(stepSize, stepSize + this.fftFrameSize_);
			for (k = 0; k < this.fftFrameSize_; k++)
				this.gOutputAccum[k] = tempArray[k];
			//Shift the input FIFO
			//These memory shifts have to be optimized.
			for (k = 0; k < inFifoLatency; k++)
				this.gInFIFO[k] = this.gInFIFO[k + stepSize];
		}
	} //for numSampsToProcess
} //process

/*
Pitchshift.prototype.processMulti = function(buffer, pitchShiftArray, osamp) {
	function setArray(array, length, val){
		var intLength = Math.floor(length);
		for (var i = 0; i < intLength; i++)
			array[i] = val;
	} //setArray
	//Prepare array of buffer clones
	var shiftedBuffers = [];
	for(var p=0; p<pitchShiftArray.length; p++)
		shiftedBuffers[p] = buffer.clone();
	//For each channel of main buffer
	for(var i=0; i<buffer.numberOfChannels; i++){
		var indata = buffer.getChannelData(i);
		numSampsToProcess = indata.length;

		//PitchShift: factor value which is between 0.5 (one octave down) and 2. (one octave up).
		var fftFrameSize2 = this.fftFrameSize_ / 2, stepSize = this.fftFrameSize_ / osamp,
			freqPerBin = this.sampleRate_ / this.fftFrameSize_, expct = 2. * Math.PI * stepSize / this.fftFrameSize_,
			inFifoLatency = this.fftFrameSize_ - stepSize, j, k = 0, magn, phase, tmp, qpd, index, signal;
		if (this.gRover_ === false)
			this.gRover_ = inFifoLatency;
		//Main processing loop
		for (j = 0; j < numSampsToProcess; j++) {
			//As long as we have not yet collected enough data just read in
			this.gInFIFO[this.gRover_] = indata[j];
			this.outdata[j] = this.gOutFIFO[this.gRover_ - inFifoLatency];
			this.gRover_++;

			//We have enough data for processing
			if (this.gRover_ >= this.fftFrameSize_) {
				this.gRover_ = inFifoLatency;
				//Windowing
				for (k = 0; k < this.fftFrameSize_; k++) {
					//Need the signal for the FFT.
					this.gFFTworksp[k] = this.gInFIFO[k] * this.hannWindow_[k];
					//this.gFFTworksp[k][1] = 0.;
				}
				this.fft.forward(this.gFFTworksp);

				//Analysis
				for (k = 0; k <= fftFrameSize2; k++) {
					//These ifs make the pitchshifter code dependent on the DFT implementation; we should decorate DFTs instead.
					if (this.algo === "FFT") {
						//Taking some "private" member out of fft here.
						magn = 2 * Math.sqrt(this.fft.real[k] * this.fft.real[k] + this.fft.imag[k] * this.fft.imag[k]);
						//aka magn = spectrum[k];
						phase = Math.atan2(this.fft.imag[k], this.fft.real[k]);
					} else if (this.algo === "RFFT") {
						//Because having the same interface but a different output schema
						//in the same library is a great fucking idea!
						var imaginary, real;
						real = this.fft.trans[k];
						imaginary = k==0 ? 0 : this.fft.trans[this.fftFrameSize_ - k];
						magn = 2 * Math.sqrt(real * real + imaginary * imaginary);
						phase = Math.atan2(imaginary, real);
					} else {
						throw new Error("Invalid DFT algorithm selected " + this.algo);
					}
					//Compute phase difference
					tmp = phase - this.gLastPhase[k];
					this.gLastPhase[k] = phase;
					tmp -= k * expct;
					qpd = tmp / Math.PI;
					if (qpd >= 0) {
						qpd = Math.floor(qpd);
						qpd += qpd & 1;
					} else {
						qpd = Math.ceil(qpd);
						qpd -= qpd & 1;
					}
					tmp -= Math.PI * qpd;
					tmp = osamp * tmp / (2 * Math.PI);
					tmp = k * freqPerBin + tmp * freqPerBin;
					this.gAnaMagn[k] = magn;
					this.gAnaFreq[k] = tmp;
				}

				for(var ppp=0; ppp < pitchShiftArray.length; ppp++){
					shiftedBuffers[ppp].getChannelData(i)[j] = shifter.outdata[j];

					//Pitch shifting
					setArray(this.gSynMagn, this.fftFrameSize_, 0);
					setArray(this.gSynFreq, this.fftFrameSize_, 0);
					for (k = 0; k <= fftFrameSize2; k++) {
						index = Math.floor(k * pitchShiftArray[ppp]);
						if (index <= fftFrameSize2) {
							this.gSynMagn[index] += this.gAnaMagn[k];
							this.gSynFreq[index] = this.gAnaFreq[k] * pitchShiftArray[ppp];
						}
					}
					//Synthesis
					for (k = 0; k <= fftFrameSize2; k++) {
						magn = this.gSynMagn[k];
						tmp = this.gSynFreq[k];
						tmp -= k * freqPerBin;
						tmp /= freqPerBin;
						tmp = 2. * Math.PI * tmp / osamp;
						tmp += k * expct;
						this.gSumPhase[k] += tmp;
						phase = this.gSumPhase[k];
						this.real_[k] = magn * Math.cos(phase);
						this.imag_[k] = magn * Math.sin(phase);
					}
					//Zero negative frequencies
					for (k = ((fftFrameSize2) + 1); (k < this.fftFrameSize_); k++) {
						this.real_[k] = 0;
						this.imag_[k] = 0;
					}
					//Do the Inverse transform
					signal = this.invFFT.inverse(this.real_, this.imag_);
					for (k = 0; k < this.fftFrameSize_; k++)
						this.gOutputAccum[k] += this.hannWindow_[k] * signal[k];
					for (k = 0; k < stepSize; k++)
						this.gOutFIFO[k] = this.gOutputAccum[k];
					var tempArray = this.gOutputAccum.slice(stepSize, stepSize + this.fftFrameSize_);
					for (k = 0; k < this.fftFrameSize_; k++)
						this.gOutputAccum[k] = tempArray[k];
					for (k = 0; k < inFifoLatency; k++)
						this.gInFIFO[k] = this.gInFIFO[k + stepSize];
				} //for pitchShiftArray
			} //if
		} //for numSampsToProcess
	} //for channels
} //processMulti
//*/
/*	//TEST
B = BUFFERS[INSTRUMENT];
var shifter = new Pitchshift( 2048, B.sampleRate, 'FFT' );
SHIFTED_BUFFERS = shifter.processMulti(B, [1/8, 1/4, 1/2, 1, 2, 4, 8], 4);
//*/

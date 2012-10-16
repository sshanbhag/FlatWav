
%
% Proposed order for calibration:
% 
% 
% 1) load xfer function
% 		
% 2) load .wav file
% 
% 3) determine max/min of xfer function


% frequencies
Fs = 500000;
Fmin = 4000;
Fmax = 80000;
stimdur = 250;
corr_frange = [4000 80000];
% hi pass filter Fc
InputHPFc = 20;
InputLPFc = 20000;
forder = 3;

SweepDuration = 300;

MIN_DB = -120;

CALLPATH = '~/Work/Data/Audio/BatCalls';
callname = 'app1.wav';


%% load xfer function data
load('test_16Oct2012_3K-110K.cal', '-MAT', 'caldata');

% plot
figure(1)
plot(caldata.freq*0.001, caldata.mag(1, :), '.-')
xlabel('Frequency (kHz)')
ylabel('dB SPL')
title('System Transfer Function')


%% synthesize test noise from min(F) to max(F)
s = synmononoise_fft(stimdur, Fs, Fmin, Fmax, 1, 0);
s = normalize(s);

% s = synmonosine(stimdur, Fs, Fsine, 1, 0);
% plot spectrum of s
figure(2)
tvec = 1000 * (0:(length(s)-1)) ./ Fs;
subplot(221)
plot(tvec, s)
title('Test signal')
xlabel('time (milliseconds)')
ylabel('V')

subplot(222)
[fraw, magraw] = daqdbfft(s, Fs, length(s));
plot(0.001*fraw, magraw);
title('Test Signal Spectrum')
xlabel('freq (kHz)')
ylabel('dB')
ylim([-120 -40])
xlim([0 0.001*Fs/2])

%% use compensate signal to compensate
[sadj, Sfull, Hnorm, foutadj] = compensate_signal(s, caldata.freq, caldata.mag(1, :), Fs, ...
							corr_frange, 'Method', 'comp', 'Normalize', 'off', 'Lowcut', 'off');

% plot compensated signal
figure(2)
subplot(223)
plot(tvec, sadj)
title('Compensated Signal')
xlabel('time (milliseconds)')
ylabel('V')

subplot(224)
[fadj, magadj] = daqdbfft(sadj, Fs, length(sadj));
plot(0.001*fadj, magadj);
title('Compensated Signal Spectrum')
xlabel('freq (kHz)')
ylabel('dB')
ylim([-120 -40])
xlim([0 0.001*Fs/2])

%% initialize NI hardware

iodev = ni_ioinit('Dev1', SweepDuration, Fs, 5);

pause(1)

%------------------------------------------------------------------------
%------------------------------------------------------------------------
%% Define a bandpass filter for processing the data
%------------------------------------------------------------------------
%------------------------------------------------------------------------
% Nyquist frequency
fnyq = iodev.Fs / 2;
% passband definition
fband = [InputHPFc InputLPFc] ./ fnyq;
% filter coefficients using a butterworth bandpass filter
[fcoeffb, fcoeffa] =	butter(forder, fband, 'bandpass');

SweepPoints = ms2samples(SweepDuration, iodev.Fs);

%% play/record raw signal
s = sin2array(s, 1, iodev.Fs);
stemp = [s; zeros(size(s))];
[resp, indx] = nidaq_calibration_io(iodev, stemp, SweepPoints);

rawresp = filter(fcoeffb, fcoeffa, sin2array(resp{1}, 5, iodev.Fs));
fftdbplot(rawresp, iodev.Fs, figure(3));

pause(1)

%% play/record raw signal
sadj = sin2array(sadj, 1, iodev.Fs);
stemp = [sadj; zeros(size(sadj))];
[resp, indx] = nidaq_calibration_io(iodev, stemp, SweepPoints);
adjresp = filter(fcoeffb, fcoeffa, sin2array(resp{1}, 5, iodev.Fs));
fftdbplot(adjresp, iodev.Fs, figure(4));


%% delete/clean  up NI subsystem

iodev = ni_ioexit(iodev); clear iodev

return





%% test with real signal

% load signal
[b, Fs, nbits, opts] = wavread(fullfile(CALLPATH, callname));
b = normalize(b)';
b = sin2array(b, 1, Fs);

% apply correction (BOOST method)
[badj, Bfull, Hnorm, foutadj] = compensate_signal(b,  caldata.freq, caldata.mag(1, :), Fs, corr_frange);

% plot call and spectrum
[fraw, magraw] = daqdbfft(b, Fs, length(b));
figure(3)
tvec = 1000 * (0:(length(b)-1)) ./ Fs;
subplot(221)
plot(tvec, b)
title('Call')
ylabel('V')

subplot(222)
plot(0.001*fraw, magraw);
title('Call Spectrum')
ylabel('dB')
ylim([-140 -40])
xlim([0 0.001*Fs/2])

% plot compensated signal
[fadj, magadj] = daqdbfft(badj, Fs, length(badj));
subplot(223)
plot(tvec, badj)
title('Compensated Call BOOST')
xlabel('time (milliseconds)')
ylabel('V')

subplot(224)
plot(0.001*fadj, magadj);
title('Compensated Call Spectrum BOOST')
xlabel('freq (kHz)')
ylabel('dB')
ylim([-140 -40])
xlim([0 0.001*Fs/2])

%% test atten method

% apply correction (BOOST method)
[aadj, Afull, Hnorm, foutadj] = compensate_signal(b,  caldata.freq, caldata.mag(1, :), Fs, corr_frange, 'ATTEN');

figure(4)
% plot call and spectrum
[fraw, magraw] = daqdbfft(b, Fs, length(b));

tvec = 1000 * (0:(length(b)-1)) ./ Fs;
subplot(221)
plot(tvec, b)
title('Call')
ylabel('V')

subplot(222)
plot(0.001*fraw, magraw);
title('Call Spectrum')
ylabel('dB')
ylim([-140 -40])
xlim([0 0.001*Fs/2])

% plot compensated signal
[fadj, magadj] = daqdbfft(aadj, Fs, length(aadj));
subplot(223)
plot(tvec, aadj)
title('Compensated Call ATTEN')
xlabel('time (milliseconds)')
ylabel('V')

subplot(224)
plot(0.001*fadj, magadj);
title('Compensated Call Spectrum ATTEN')
xlabel('freq (kHz)')
ylabel('dB')
ylim([-140 -40])
xlim([0 0.001*Fs/2])


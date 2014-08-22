%------------------------------------------------------------------------
%------------------------------------------------------------------------
% script to exercise the compensate_signal.m function
%
% Proposed order for calibration:
% 
% 
% 1) load xfer function
% 		
% 2) load .wav file
% 
% 3) determine max/min of xfer function
%------------------------------------------------------------------------
%------------------------------------------------------------------------

%------------------------------------------------------------------------
% CONSTANTS
%------------------------------------------------------------------------
% minimum dB level
MIN_DB = -120;

%------------------------------------------------------------------------
%% settings
%------------------------------------------------------------------------
% sample rate
Fs = 500000;
Fmin = 100;
Fmax = 100000;
% signal duration
stimdur = 250;
% [lower  upper] frequency bounds for compensation
corr_frange = [5000 100000];
% low freq cutoff level
lowcut_level = 4000;
% normalization value
norm_level = 1;
% db adj level (for compress compensation method)
dbadj_level = 65;

% path to calls, call name
CALLPATH = '~/Work/Data/Audio/BatCalls/BatStrings_Emily';
callname = 'app1.wav'
% path to calibration data, calibration file
CALPATH = '~/Work/Data/Audio/Calibration/MouseRig_21August2014';
calfile = '21Aug2014_working.cal'


%------------------------------------------------------------------------
%% load xfer function data
%------------------------------------------------------------------------
load(fullfile(CALPATH, calfile), '-MAT', 'caldata');

% smoothed data
smoothed = smooth_calibration_data(1, caldata, 5);
% plot
figure(1)
plot(	0.001*caldata.freq, caldata.mag(1, :), 'k.-', ...
		0.001*caldata.freq, smoothed(1, :), 'b.-')
xlabel('Frequency (kHz)')
ylabel('dB SPL')
title('System Transfer Function')

%------------------------------------------------------------------------
%% synthesize test noise from min(F) to max(F)
%------------------------------------------------------------------------
s = synmononoise_fft(stimdur, Fs, Fmin, Fmax, 1, 0);
s = normalize(s);

% s = synmonosine(stimdur, Fs, Fsine, 1, 0);
% plot spectrum of s
[fraw, magraw] = daqdbfft(s, Fs, length(s));
figure(2)
tvec = 1000 * (0:(length(s)-1)) ./ Fs;
subplot(221)
plot(tvec, s)
title('Test signal')
xlabel('time (milliseconds)')
ylabel('V')

subplot(222)
plot(0.001*fraw, magraw);
title('Test Signal Spectrum')
xlabel('freq (kHz)')
ylabel('dB')
ylim([-120 -40])
xlim([0 0.001*Fs/2])

%% use compensate signal to compensate
[sadj, Sfull, Hnorm, foutadj] = ...
						compensate_signal( ...
									s, ...
									caldata.freq, caldata.mag(1, :), ...
									Fs, ...
									corr_frange, ...
									'Method', 'boost', ...
									'Lowcut', lowcut_level, ...
 									'Normalize', 'off', ... 
									'Level', dbadj_level);

% plot compensated signal
[fadj, magadj] = daqdbfft(sadj, Fs, length(sadj));
figure(2)
subplot(223)
plot(tvec, sadj)
title('Compensated Signal')
xlabel('time (milliseconds)')
ylabel('V')

subplot(224)
plot(0.001*fadj, magadj);
title('Compensated Signal Spectrum')
xlabel('freq (kHz)')
ylabel('dB')
ylim([-120 -40])
xlim([0 0.001*Fs/2])



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


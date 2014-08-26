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
Fs = 250000;
Fmin = 100;
Fmax = 100000;
% signal duration
stimdur = 250;
% [lower  upper] frequency bounds for compensation
corr_frange = [4000 100000];
% low freq cutoff level
lowcut_level = 3000;
% normalization value
norm_level = 1;
% db adj level (for compress compensation method)
dbadj_level = 65;
% compensation method
compmethod = 'compress';
% pre-filter the signal before compensation?
prefilter = 'on';
% filter the signal after compensation?
% postfilter = [5000 95000];
postfilter = 'off';
% limit calculation of compensation values to the range being corrected?
range_limit = 'on';
% put a cap on the amount of correction?
corr_limit = 'off';
% smooth edges at boundaries of corr_frange before performing iFFT?
smooth_edges = 'off';

specwin = 512;

% path to calls, call name
CALLPATH = pwd;
callname = 'app1.wav'
% path to calibration data, calibration file
CALPATH = pwd;
calfile = '22Aug2014_test_CapOn.cal'


%------------------------------------------------------------------------
%% load xfer function data
%------------------------------------------------------------------------
load(fullfile(CALPATH, calfile), '-MAT', 'caldata');

% smoothed data
smoothed = smooth_calibration_data(1, caldata, 3);

%{
% plot
figure(1)
plot(	0.001*caldata.freq, caldata.mag(1, :), 'k.-', ...
		0.001*caldata.freq, smoothed(1, :), 'b.-')
xlabel('Frequency (kHz)')
ylabel('dB SPL')
title('System Transfer Function')
%}

%{
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
									caldata.freq, smoothed(1, :), ...
									Fs, ...
									corr_frange, ...
									'Method', compmethod, ...
									'Lowcut', lowcut_level, ...
 									'Normalize', norm_level, ... 
									'Level', dbadj_level, ...
									'Prefilter', prefilter);

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
title(sprintf('Compensated Signal Spectrum %s', compmethod))
xlabel('freq (kHz)')
ylabel('dB')
ylim([-120 -40])
xlim([0 0.001*Fs/2])

%% spectrogram plot
figure(3)

subplot(211)
[S, F, T, P] = spectrogram(	s, ...
										specwin, ...
										[], ...
										specwin, ...
										Fs	);
P = 20*log10(P);
P(P == -Inf) = min(min(P(P ~= -Inf)));	
surf(1000*T, 0.001*F, P, 'edgecolor', 'none');
view(0, 90);
xlim([min(tvec) max(tvec)])
ylim([0 0.001*Fs/2]);
colormap(gray)

subplot(212)
[S, F, T, P] = spectrogram(	sadj, ...
										specwin, ...
										[], ...
										specwin, ...
										Fs	);
P = 20*log10(P);
P(P == -Inf) = min(min(P(P ~= -Inf)));	
surf(1000*T, 0.001*F, P, 'edgecolor', 'none');
view(0, 90);
xlim([min(tvec) max(tvec)])
ylim([0 0.001*Fs/2]);
xlabel('Time (ms)')
colormap(gray)
% 	caxis([min(min(P)) max(max(P))])
%}

%------------------------------------------------------------------------
%% test with real signal
%------------------------------------------------------------------------

% load signal
[b, Fs, nbits, opts] = wavread(fullfile(CALLPATH, callname));
% b = normalize(b)';
b = b';
b = sin2array(b, 1, Fs);

% apply correction (BOOST method)
[badj, Bfull, Hnorm, foutadj] = compensate_signal( ...
											b, ...
											caldata.freq, smoothed(1, :), ...
											Fs, ...
											corr_frange, ...
											'Method', compmethod, ...
											'Lowcut', lowcut_level, ...
											'Normalize', norm_level, ... 
											'Level', dbadj_level, ...
											'Postfilter', postfilter, ...
											'Prefilter', prefilter, ...
											'Rangelimit', range_limit, ...
											'Corrlimit', corr_limit, ...
											'SmoothEdges', smooth_edges	);
% plot call and spectrum
[fraw, magraw] = daqdbfft(b, Fs, length(b));
figure
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
title(sprintf('Compensated Call Spectrum %s', compmethod))
xlabel('freq (kHz)')
ylabel('dB')
ylim([-140 -40])
xlim([0 0.001*Fs/2])

% spectrogram plot
figure
subplot(211)
[S, F, T, P] = spectrogram(	b, ...
										specwin, ...
										[], ...
										specwin, ...
										Fs	);
P = 20*log10(P);
P(P == -Inf) = min(min(P(P ~= -Inf)));	
surf(1000*T, 0.001*F, P, 'edgecolor', 'none');
view(0, 90);
xlim([min(tvec) max(tvec)])
ylim([0 0.001*Fs/2]);
colormap(gray)

subplot(212)
[S, F, T, P] = spectrogram(	badj, ...
										specwin, ...
										[], ...
										specwin, ...
										Fs	);
P = 20*log10(P);
P(P == -Inf) = min(min(P(P ~= -Inf)));	
surf(1000*T, 0.001*F, P, 'edgecolor', 'none');
view(0, 90);
xlim([min(tvec) max(tvec)])
ylim([0 0.001*Fs/2]);
xlabel('Time (ms)')
colormap(gray)
% 	caxis([min(min(P)) max(max(P))])


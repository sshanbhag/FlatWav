
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
Fmin = 100;
Fmax = 100000;
stimdur = 250;
corr_frange = [5000 100000];

MIN_DB = -120;

CALLPATH = '~/Work/Data/Audio/BatCalls';
callname = 'app1.wav';


%% load xfer function data
load('TDT3972_5V_MicNoScreen_cal.mat', '-MAT', 'caldata');

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
[sadj, Sfull, Hnorm, foutadj] = compensate_signal(s, caldata, Fs, corr_frange);

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



%% test signal


% load signal
[b, Fs, nbits, opts] = wavread(fullfile(CALLPATH, callname));
b = normalize(b)';
b = sin2array(b, 1, Fs);

% apply correction (BOOST method)
[badj, Bfull, Hnorm, foutadj] = compensate_signal(b, caldata, Fs, corr_frange);

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
title('Compensated Call')
xlabel('time (milliseconds)')
ylabel('V')

subplot(224)
plot(0.001*fadj, magadj);
title('Compensated Call Spectrum')
xlabel('freq (kHz)')
ylabel('dB')
ylim([-140 -40])
xlim([0 0.001*Fs/2])


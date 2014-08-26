%------------------------------------------------------------------------
%------------------------------------------------------------------------
% script to test fft/ifft synverse method
%------------------------------------------------------------------------
%------------------------------------------------------------------------

%------------------------------------------------------------------------
% CONSTANTS
%------------------------------------------------------------------------
% arbitrary minimum dB value
MIN_DB = -120;
% need to have a small, but non-zero value when taking log, so set that here
ZERO_VAL = 1e-17;

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
corr_frange = [5000 100000];
% low freq cutoff level
lowcut_level = 5000;
% normalization value
norm_level = 1;
% db adj level (for compress compensation method)
dbadj_level = 200;
specwin = 512;
compmethod = 'boost';
prefilter = 'on';

% path to calls, call name
CALLPATH = pwd;
callname = 'app1.wav';


%------------------------------------------------------------------------
%% synthesize test noise from min(F) to max(F)
%------------------------------------------------------------------------
s = synmononoise_fft(stimdur, Fs, Fmin, Fmax, 1, 0);
s = normalize(s);

% get spectrum of s signal
% length of signal
slen = length(s);
% for speed's sake, get the nearest power of 2 to the desired output length
NFFT = 2.^(nextpow2(slen));
% fft
S = fft(s, NFFT);
%non-redundant points are kept
Nunique = NFFT/2;
Sunique = S(1:Nunique);
% get the magnitudes of the FFT  and scale by 2 because we're taking only
% half of the points from the "full" FFT vector S;
Smag = abs(Sunique)/slen;
Smag(2:end) = 2*Smag(2:end);
% get phase
Sphase = angle(Sunique);
% convert to db - need to avoid log(0)
tmp = Smag;
tmp(tmp==0) = ZERO_VAL; 
SdBmag = db(tmp);
clear tmp;

% convert back to linear scale...
Sadj = invdb(SdBmag);

% scale for length of signal and divide by 2 to scale for conversion to 
% full FFT before inverse FFT
Sadj = slen * Sadj;
%Sadj = slen * Sadj ./ 2;
%Sadj = NFFT * Sadj ./ 2;
%Sadj = NFFT * Sadj;

% create compensated time domain signal from spectrum
[sadj, Sfull] = synverse(Sadj, Sphase, 'DC', 'no');
% return only 1:slen points
sadj = sadj(1:slen);

% plot s and spectrum of s
[fraw, magraw] = daqdbfft(s, Fs, length(s));
figure(1)
tvec = 1000 * (0:(length(s)-1)) ./ Fs;
subplot(231)
plot(tvec, s)
title('Test signal')
xlabel('time (milliseconds)')
ylabel('V')
ylim(max(abs(s)) * [-1 1])
subplot(232)
plot(0.001*fraw, magraw);
title('Test Signal Spectrum')
xlabel('freq (kHz)')
ylabel('dB')
ylim([-120 -40])
xlim([0 0.001*Fs/2])

% spectrogram plot
subplot(233)
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

% plot s, sadj 
[fadj, magadj] = daqdbfft(sadj, Fs, length(sadj));
subplot(234)
plot(tvec, sadj)
ylim(max(abs(s)) * [-1 1])
title('Compensated Signal')
xlabel('time (milliseconds)')
ylabel('V')

subplot(235)
plot(0.001*fadj, magadj);
title('processed test signal')
xlabel('freq (kHz)')
ylabel('dB')
ylim([-120 -40])
xlim([0 0.001*Fs/2])

subplot(236)
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

%------------------------------------------------------------------------
%% test with real signal
%------------------------------------------------------------------------

% load signal
[b, Fs, nbits, opts] = wavread(fullfile(CALLPATH, callname));
% b = normalize(b)';
b = b';
b = sin2array(b, 0.5, Fs);

% get spectrum of s signal
% length of signal
blen = length(b);
% for speed's sake, get the nearest power of 2 to the desired output length
NFFT = 2.^(nextpow2(blen));
% fft
B = fft(b, NFFT);
%non-redundant points are kept
Nunique = NFFT/2;
Bunique = B(1:Nunique);
% get the magnitudes of the FFT  and scale by 2 because we're taking only
% half of the points from the "full" FFT vector S;
Bmag = abs(Bunique)/blen;
Bmag(2:end) = 2*Bmag(2:end);
% get phase
Bphase = angle(Bunique);
% convert to db - need to avoid log(0)
tmp = Bmag;
tmp(tmp==0) = ZERO_VAL; 
BdBmag = db(tmp);
clear tmp;

% convert back to linear scale...
Badj = invdb(BdBmag);

% scale for length of signal and divide by 2 to scale for conversion to 
% full FFT before inverse FFT
Badj = blen * Badj;
%Badj = blen * Badj ./ 2;
%Badj = NFFT * Badj ./ 2;
%Badj = NFFT * Badj;

% create compensated time domain signal from spectrum
[badj, Bfull] = synverse(Badj, Bphase, 'DC', 'no');
% return only 1:blen points
badj = badj(1:blen);

% plot b and spectrum of b
[fraw, magraw] = daqdbfft(b, Fs, length(b));
figure(2)

% raw signal plot
tvec = 1000 * (0:(length(b)-1)) ./ Fs;
subplot(231)
plot(tvec, b)
title('WAV signal')
xlabel('time (milliseconds)')
ylabel('V')
ylim(max(abs(b)) * [-1 1])
% raw signal magnitude spectrum plot
subplot(232)
plot(0.001*fraw, magraw);
title('WAV Spectrum')
xlabel('freq (kHz)')
ylabel('dB')
ylim([-120 -40])
xlim([0 0.001*Fs/2])
% raw signal spectrogram plot
subplot(233)
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


% plot b, badj 
[fadj, magadj] = daqdbfft(badj, Fs, length(badj));
% adj signal plot
subplot(234)
plot(tvec, badj)
ylim(max(abs(b)) * [-1 1])
title('Compensated Signal')
xlabel('time (milliseconds)')
ylabel('V')
% adj mag spectrum plot
subplot(235)
plot(0.001*fadj, magadj);
title('processed test signal')
xlabel('freq (kHz)')
ylabel('dB')
ylim([-120 -40])
xlim([0 0.001*Fs/2])
% adj spectrogram plot
subplot(236)
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


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
Fmin = 5000;
Fmax = 115000;
stimdur = 500;
Fsine = 10000;
corr_frange = [5000 110000];

MIN_DB = -120;

CALLPATH = '~/Work/Data/Audio/BatCalls';
callname = 'app1.wav';

%% synthesize test noise from min(F) to max(F)
s = synmononoise_fft(stimdur, Fs, Fmin, Fmax, 1, 0);
s = normalize(s);
% s = synmonosine(stimdur, Fs, Fsine, 1, 0);
tvec = 1000 * (0:(length(s)-1)) ./ Fs;
figure(2)
subplot(211)
plot(tvec, s)
title('test signal')
xlabel('time (milliseconds)')
ylabel('V')

% get fft of test noise

% length of signal
Nsignal = length(s);
% for speed's sake, get the nearest power of 2 to the desired output length
NFFT = 2.^(nextpow2(Nsignal));
% fft
S = fft(s, NFFT);
%non-redundant points are kept
Nunique = NFFT/2;
Sunique = S(1:Nunique);
% get the magnitudes of the FFT  and scale by 2 because we're taking only
% half of the points from the "full" FFT vector S;
Smag = abs(Sunique)/Nsignal;
Smag(2:end) = 2*Smag(2:end);
% get phase
Sphase = angle(Sunique);
% convert to db (need to avoid log(0))
tmp = Smag;
tmp(tmp==0) = 1e-17; 
SdBmag = db(tmp);
clear tmp;
% build frequency vector
% This is an evenly spaced frequency vector with Nunique points.
% scaled by the Nyquist frequency (Fn ==1/2 sample freq.)
f = (Fs/2)*linspace(0, 1, NFFT/2);
% plot spectrum
figure(2)
subplot(212)
plot(0.001*f, SdBmag);
title('Test Signal Spectrum')
xlabel('freq (kHz)')
ylabel('dB')


%% load xfer function data
load('TDT3972_5V_MicNoScreen_cal.mat', '-MAT', 'caldata');

% pull out freq and magnitude (SPL) data
F = caldata.freq;
H = caldata.mag(1, :);
% plot
figure(1)
subplot(211)
plot(F*0.001, H, '.-')
xlabel('Frequency (kHz)')
ylabel('dB SPL')
title('xfer function')

% normalize to peak of xfer function
% find peak, peak index
[hpeak, hpeak_ind] = max(H)
% get freq at peak
fpeak = F(hpeak_ind);
% normalize by finding deviation from peak
Hnorm = hpeak - H;
% plot
subplot(212)
plot(F*0.001, Hnorm, '.-')
xlabel('Frequency (kHz)')
ylabel('normalized dB SPL')
title('Normalized Xfer')

%% apply correction (BOOST method)
% 
% procedure:	find compensation values for frequency range for which there is
%					calibration data and apply to FFT, then iFFT to get corrected version

if corr_frange(1) < min(caldata.freq)
	corr_frange(1) = min(caldata.freq);
end
if corr_frange(2) > max(caldata.freq)
	corr_frange(2) = max(caldata.freq);
end

% first, need to find max, min of calibration range
valid_indices = find(between(f, corr_frange(1), corr_frange(2))==1);
% check to make sure there is overlap in ranges
if isempty(valid_indices)
	error('%s: mismatch between FFT frequencies and calibration range', mfilename);
end
% then, get the frequencies for that range
corr_f = f(valid_indices);
% interpolate to get the correction values (in dB!)
corr_vals = interp1(F, Hnorm, corr_f);

% create adjusted magnitude vector from Smag (in dB)
SdBadj = SdBmag;
% apply correction
SdBadj(valid_indices) = SdBadj(valid_indices) + corr_vals;

% set freqs below 4KHz to MINDB
sub4indices = find(f < 4000);
if ~isempty(sub4indices)
	SdBadj(sub4indices) = MIN_DB * ones(size(SdBadj(sub4indices)));
end

% plot spectrum
figure(4)
subplot(211)
plot(0.001*f, SdBadj)
title('Compensated Signal Spectrum')
xlabel('freq (kHz)')
ylabel('dB')

% convert back to linear scale...
Sadj = invdb(SdBadj);
% and plot
subplot(212)
plot(0.001*f, Sadj)
title('Compensated Signal Spectrum')
xlabel('freq (kHz)')
ylabel('arbitrary values')

% scale for length of signal and divide by 2 to scale for conversion to 
% full FFT before inverse FFT
Sadj = Nsignal * Sadj ./ 2;
% create compensated time domain signal from spectrum
[sadj, Sfull] = synverse(Sadj, Sphase, 'DC', 'no');
sadj = normalize(sadj(1:Nsignal));

% plot compensated signal
[f2, mag2] = daqdbfft(sadj, Fs, length(sadj));
figure(5)
subplot(311)
plot(s)
subplot(312)
plot(sadj);
subplot(313)
plot(0.001*f2(2:end), mag2(2:end))



%% test signal


% load signal
[b, Fs, nbits, opts] = wavread(fullfile(CALLPATH, callname));
b = b';
b = sin2array(b, 1, Fs);
% get fft of test signal
% length of signal
Nsignal = length(b);
% for speed's sake, get the nearest power of 2 to the desired output length
NFFT = 2.^(nextpow2(Nsignal));
% fft
B = fft(b, NFFT);
%non-redundant points are kept
Nunique = NFFT/2;
Bunique = B(1:Nunique);
% get the magnitudes of the FFT  and scale by 2 because we're taking only
% half of the points from the "full" FFT vector S;
Bmag = abs(Bunique)/Nsignal;
Bmag(2:end) = 2*Bmag(2:end);
% get phase
Bphase = angle(Bunique);
% convert to db (need to avoid log(0))
tmp = Bmag;
tmp(tmp==0) = 1e-17; 
BdBmag = db(tmp);
clear tmp;
% build frequency vector
% This is an evenly spaced frequency vector with Nunique points.
% scaled by the Nyquist frequency (Fn ==1/2 sample freq.)
f = (Fs/2)*linspace(0, 1, NFFT/2);
% plot spectrum
figure(6)
subplot(221)
plot(b);
subplot(223)
plot(0.001*f, BdBmag);
title('Call Spectrum')
xlabel('freq (kHz)')
ylabel('dB')


% apply correction (BOOST method)
% first, need to find max, min of calibration range
valid_indices = find(between(f, corr_frange(1), corr_frange(2))==1);
% check to make sure there is overlap in ranges
if isempty(valid_indices)
	error('%s: mismatch between FFT frequencies and calibration range', mfilename);
end
% then, get the frequencies for that range
corr_f = f(valid_indices);
% interpolate to get the correction values (in dB!)
corr_vals = interp1(F, Hnorm, corr_f);

% create adjusted magnitude vector from Smag (in dB)
BdBadj = BdBmag;
% apply correction
BdBadj(valid_indices) = BdBadj(valid_indices) + corr_vals;

% set freqs below 4KHz to MINDB
sub4indices = find(f < 4000);
if ~isempty(sub4indices)
	BdBadj(sub4indices) = MIN_DB * ones(size(BdBadj(sub4indices)));
end

% convert back to linear scale...
Badj = invdb(BdBadj);
% scale for length of signal and divide by 2 to scale for conversion to 
% full FFT before inverse FFT
Badj = Nsignal * Badj ./ 2;
% create compensated time domain signal from spectrum
[badj, Bfull] = synverse(Badj, Bphase, 'DC', 'no');
badj = normalize(badj(1:Nsignal));
% plot compensated signal
[f2, bmag2] = daqdbfft(badj, Fs, length(badj));

figure(6)
subplot(222)
plot(badj)
title('compensated signal')
subplot(224)
plot(0.001*f2, bmag2);
title('Compensated Spectrum')
xlabel('freq (kHz)')




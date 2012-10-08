% frequencies
Fs = 44100;
Fmin = 100;
Fmax = 10000;
stimdur = 500;
corr_frange = [200 10000];

%------------------------------------------------------------------------
%------------------------------------------------------------------------
% define some constants
%------------------------------------------------------------------------
%------------------------------------------------------------------------
MIN_DB = -120;
ZERO_VAL = 1e-17;


%% load dummy calibration  data
caldata = fake_caldata('FREQS', 100:100:10000);

% assign meaningful values to caldata.mag
x = 2*pi*linspace(0, 1, length(caldata.mag(1, :)));
m = 20 * sin(x) + 40;
caldata.mag = [m; m];

% shorter names for freq and mag from caldata
F = caldata.freq;
M = caldata.mag(1, :);

% plot
figure(1)
subplot(211)
plot(0.001*F, M, '.-')
xlabel('Frequency (kHz)')
ylabel('dB SPL')
title('System Transfer Function')

%% synthesize test noise from min(F) to max(F)
s = synmononoise_fft(stimdur, Fs, Fmin, Fmax, 1, 0);
s = normalize(s);
s = sin2array(s, 5, Fs);

% s = synmonosine(stimdur, Fs, Fsine, 1, 0);
% plot spectrum of s
[fraw, magraw, phiraw] = daqdbfullfft(s, Fs, length(s));
figure(2)
tvec = 1000 * (0:(length(s)-1)) ./ Fs;
subplot(231)
plot(tvec, s)
title('Test signal')
xlabel('time (milliseconds)')
ylabel('V')

subplot(232)
plot(0.001*fraw, magraw);
title('Test Signal Magnitude')
xlabel('freq (kHz)')
ylabel('dB')
ylim([-120 0])
xlim([0 0.001*Fs/2])

subplot(233)
plot(0.001*fraw, phiraw);
title('Test Signal Phase');
xlabel('freq (kHz)');
ylabel('rad');
xlim([0 0.001*Fs/2]);

% 
% 
% % play sound
% sound(0.9*s, Fs);
% pause(1);


%% test compress method
%[sadj, Sfull, Hnorm, foutadj] = compensate_signal(s, F, M, Fs, corr_frange, 'Method', 'BOOST', 'Normalize', 'on', 'Lowcut', 'off');

calfreq = F;
calmag = M;
LOWCUT = 0;

%------------------------------------------------------------------------
%------------------------------------------------------------------------
% get spectrum of signal, s
%------------------------------------------------------------------------
%------------------------------------------------------------------------
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
% convert to db - need to avoid log(0)
tmp = Smag;
tmp(tmp==0) = ZERO_VAL; 
SdBmag = db(tmp);
clear tmp;

%------------------------------------------------------------------------
%------------------------------------------------------------------------
% determine frequency vector for calibration of signal
%------------------------------------------------------------------------
%------------------------------------------------------------------------
% build frequency vector
% This is an evenly spaced frequency vector with Nunique points.
% scaled by the Nyquist frequency (Fn ==1/2 sample freq.)
f = (Fs/2)*linspace(0, 1, NFFT/2);

% check the correction frequency range, adjust values if out of bounds
if corr_frange(1) < min(calfreq)
	corr_frange(1) = min(calfreq);
end
if corr_frange(2) > max(calfreq)
	corr_frange(2) = max(calfreq);
end

% need to find max, min of calibration range
valid_indices = find(between(f, corr_frange(1), corr_frange(2))==1);

% check to make sure there is overlap in ranges
if isempty(valid_indices)
	% if not, throw an error
	error('%s: mismatch between FFT frequencies and calibration range', mfilename);
end

% then, get the frequencies for correcting that range
corr_f = f(valid_indices);

%------------------------------------------------------------------------
%------------------------------------------------------------------------
% apply correction using COMPRESS method
%------------------------------------------------------------------------
%------------------------------------------------------------------------
% procedure:	find subtractive compensation values for frequency range 
%					for which there are calibration data and apply to FFT, 
%					then iFFT to get corrected version
%------------------------------------------------------------------------
% some assumptions:
% 					magnitude values are in ACTUAL, dB SPL range.  
% 						¡this algorithm blows up for negative magnitudes!
%------------------------------------------------------------------------
%------------------------------------------------------------------------
% find max and min in magnitude spectrum
maxmag = max(calmag);
minmag = min(calmag);
% compute middle value
midmag = ((maxmag - minmag) / 2) + minmag
meanmag = mean(calmag)
% normalize by finding deviation from peak
Magnorm = midmag - calmag;

figure(1)
subplot(212)
plot(0.001*calfreq, Magnorm, '.-')
xlabel('Frequency (kHz)')
ylabel('dB')
title('Correction Function: COMPRESS method')

% interpolate to get the correction values (in dB!)
corr_vals = interp1(calfreq, Magnorm, corr_f);

% create adjusted magnitude vector from Smag (in dB)
SdBadj = SdBmag;
% apply correction
SdBadj(valid_indices) = SdBadj(valid_indices) + corr_vals;

% set freqs below LOWCUT to MINDB
if LOWCUT > 0
	lowcutindices = find(f < LOWCUT);
	if ~isempty(lowcutindices)
		SdBadj(lowcutindices) = MIN_DB * ones(size(SdBadj(lowcutindices)));
	end
end

% convert back to linear scale...
Sadj = invdb(SdBadj);

% scale for length of signal and divide by 2 to scale for conversion to 
% full FFT before inverse FFT
Sadj = Nsignal * Sadj ./ 2;

% create compensated time domain signal from spectrum
[sadj, Sfull] = synverse(Sadj, Sphase, 'DC', 'no');
% return only 1:Nsignal points
sadj = sadj(1:Nsignal);


% plot compensated signal
[fadj, magadj, phiadj] = daqdbfullfft(sadj, Fs, length(sadj));
figure(2)
subplot(234)
plot(tvec, sadj)
title('Boost Comp Signal')
xlabel('time (milliseconds)')
ylabel('V')

subplot(235)
plot(0.001*fadj, magadj);
title('Boost Comp Magn.')
xlabel('freq (kHz)')
ylabel('dB')
ylim([-120 0])
xlim([0 0.001*Fs/2])

subplot(236)
plot(0.001*fadj, phiadj);
title('Boost Comp Phase');
xlabel('freq (kHz)');
ylabel('rad');
xlim([0 0.001*Fs/2]);


% % play sound
% sound(sin2array(0.9*sadj, 1, Fs), Fs);
% pause(1);
%}

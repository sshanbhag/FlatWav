
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
Fs = 250000;
Fmin = 1000;
Fmax = 120000;
stimdur = 100;
Fsine = 10000;

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


%% synthesize test noise from min(F) to max(F)
% s = synmononoise_fft(stimdur, Fs, Fmin, Fmax, 1, 0);
% s = normalize(s);
s = synmonosine(stimdur, Fs, Fsine, 1, 0);
tvec = 1000 * (0:(length(s)-1)) ./ Fs;
figure(2)
plot(tvec, s)
title('test signal')
xlabel('time (milliseconds)')
ylabel('V')

%% get fft of test noise

% length of signal
N = length(s);

% for speed's sake, get the nearest power of 2 to the desired output length
NFFT = 2.^(nextpow2(N));

% fft
S = fft(s, NFFT);
%non-redundant points are kept
Nunique = NFFT/2;
Sunique = S(1:Nunique);
% get the magnitudes of the FFT  and scale by 2 because we're taking only
% half of the points from the "full" FFT vector S;
Smag = abs(Sunique)/N;
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
figure(3)
plot(0.001*f, SdBmag);
title('Test Signal Spectrum')
xlabel('freq (kHz)')
ylabel('dB')


%% apply correction (BOOST method)
% 
% procedure:	find compensation values for frequency range for which there is
%					calibration data and apply to FFT, then iFFT to get corrected version

% first, need to find max, min of calibration range
corr_frange = [5000 100000];

% valid_indices = find(between(f, corr_frange(1), corr_frange(2)));
[min_diff, valid_indices] = min(abs(f - Fsine));

corr_f = f(valid_indices);
corr_vals = interp1(F, Hnorm, corr_f);

SdBadj = SdBmag;
SdBadj(valid_indices) = SdBadj(valid_indices) + corr_vals';
figure(4)
plot(0.001*f, SdBadj)
title('Compensated Signal Spectrum')
xlabel('freq (kHz)')
ylabel('dB')

% convert back to linear scale
figure(5)
% Sadj = invdb(SdBadj);
Sadj = invdb(SdBmag);
Sadj = N * Sadj ./ 2;
plot(0.001*f, Sadj)
title('Compensated Signal Spectrum')
xlabel('freq (kHz)')
ylabel('v')


% build the full, complex FFT arrays
Sred = complex(Sadj.*cos(Sphase), Sadj.*sin(Sphase));
% Sfull = buildfft(Sred);

fftred = Sred;

% N is total number of points in the spectrum
m = length(fftred);

% allocate the net spectrum fftfull
fftfull = zeros(1, 2*m);

% % first portion of fftfull is same as fftred
% % leave out the DC component (fftred(1))
% fftfull(2:m) = fftred(2:m);
% 
% % second portion is complex conjugate of Sreduced and in reverse order
% % (setting  DC component to zero which is at fftreduced(1) and fftfull(end))
% 
% fftfull((m+1):((2*m)-1)) = conj(fftred(m:-1:2));

a1 = zeros(1, m);

a1(2:m) = fftred(2:m);
a2 = fliplr(a1);
a2 = real(a2) - 1i*imag(a2);
fftfull = [a1 a2];

figure(6)

plot([abs(fftfull) - abs(S)], '.');



S1 = (S(1:Nunique));
S2 = (S((Nunique+1):end));
return





% inverse fft
sadj = ifft(fftfull);




NFFT = length(S) / 2;

index1 = 1:NFFT;
index2 = (2*NFFT):-1:(NFFT+1);

out = S(index1) .* S(index2);



% 
% 
% corr_f = F(1):stim_fft_step:F(end);
% 
% corr_mag = interp1(F, Hnorm, corr_f);





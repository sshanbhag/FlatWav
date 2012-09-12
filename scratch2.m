
% frequencies
Fs = 1000;
Fmin = 1000;
Fmax = 2000;
Fsine = 100;
duration = 1000;

%% synthesize test noise from min(F) to max(F)
% s = synmononoise_fft(100, Fs, Fmin, Fmax, 1, 0);
% s = normalize(s);
s = synmonosine(duration, Fs, Fsine, 1, 0);
% add DC for fun
s = s + 1;
tvec = 1000 * (0:(length(s)-1)) ./ Fs;


%% get fft of test noise
N = length(s);
NFFT = 2.^(nextpow2(N));
S = fft(s, NFFT);

%non-redundant points are kept
Nunique = (NFFT/2) + 1;
Sunique = S(1:Nunique);
% get the magnitudes of the FFT  and scale by 2 because we're taking only
% half of the points from the "full" FFT vector S;
Smag = abs(Sunique)/N;
Smag(2:end) = 2*Smag(2:end);
% get phase
Sphase = angle(Sunique);
% build frequency vector
% This is an evenly spaced frequency vector with Nunique points.
% scaled by the Nyquist frequency (Fn ==1/2 sample freq.)
f = (Fs/2)*linspace(0, 1, (NFFT/2)+1);

figure(2)
plot(f, abs(Sunique), '.-')
grid
title('Spectrum - Magnitude')
xlabel('Frequency (Hz)')
ylabel('Magnitude')

% now, rebuild spectrum
% properly scale the magnitude
Sadj = N * [Smag(1) Smag(2:end)/2];
% build complex spectrum
Scomplex = complex(Sadj.*cos(Sphase), Sadj.*sin(Sphase));
% assign to Snew
indx1 = 1:Nunique;
indx2 = (Nunique+1):NFFT;
Snew = zeros(size(S));
Snew(indx1) = Scomplex;
Snew(indx2) = conj(fliplr(Scomplex(2:(end-1))));

% plot the original magnitude and reconstructed magnitude
figure(3)
plot(abs(S), 'r.-')
hold on
	plot(abs(Snew), 'g.')
hold off

snew = ifft(Snew);



%% test inverse FFT
s2 = ifft(S);

figure(1)
subplot(311)
plot(tvec, s)
title('test signal')
xlabel('time (milliseconds)')
ylabel('V')

subplot(312)
plot(tvec, s2(1:N))

subplot(313)
plot(tvec, s - s2(1:N))

%% test synverse

[g, G] = synverse(Smag, Sphase, 'DC', 'yes');
figure(4)
plot(abs(G(1:10)) - abs(Snew(1:10)))

%% synverse_test.m tests synverse.m algorithm

clear all


%% test settings
Fs = 10000;
Fmin = 100;
Fmax = 3000;
Fsine = 100;
duration = 1000;
DCval = 0;

%synthesize test noise and get FFT
s = synmononoise_fft(1000, Fs, Fmin, Fmax, 1, 0);
s = normalize(s);
% s = synmonosine(duration, Fs, Fsine, 1, 0);
% add DC for fun
s = s + DCval;
% build time vector for plots
tvec = 1000 * (0:(length(s)-1)) ./ Fs;

% get fft of test noise
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

figure(1)
subplot(311)
plot(f, Smag, '.-')
grid
title('Smag')
xlabel('Frequency (Hz)')
ylabel('Magnitude')

subplot(312)
plot(abs(S), '.-')
grid
title('abs(S)')

%------------------------------------------------------------------------------
%------------------------------------------------------------------------------
%------------------------------------------------------------------------------
% synverse algorithm here
%------------------------------------------------------------------------------
%------------------------------------------------------------------------------

%% set defaults
HAS_DC = 1;
% Smag = Smag(2:end);
% Sphase = Sphase(2:end);

%% check inputs
%{
nvararg = length(varargin);
if nvararg
	aindex = 1;
	while aindex < nvararg
		switch(upper(varargin{aindex}))
			case 'DC'
				if strcmpi(varargin{aindex+1}, 'no')
					HAS_DC = 0;
				elseif strcmp(varargin{aindex+1}, 'yes')
					HAS_DC = 1;
				else
					error('%s: unknown setting %s for DC option', ...
									mfilename, varargin{aindex+1});
				end
				aindex = aindex + 2;
			otherwise
				error('%s: unknown argument %s', mfilename, varargin{aindex});
		end
	end
end
%}


%% build FFT

% scale the magnitudes for conversion to full length spectrum
if HAS_DC
	% if HAS_DC is set, the Smag vector has the DC value of the 
	% signal at Smag(1).  
	Nsignal = length(Smag) - 1;
	Nsignal = N;
	Sadj = Nsignal * [Smag(1) Smag(2:end)./2];
else
	Nsignal = length(Smag);
	Nsignal = N;
	Sadj = Nsignal * [0 Smag./2];
	Sphase = [0 Sphase];
end

% Nunique is # of unique (non-conjugate) points in the full spectrum, 
% including the DC component (or, length of the Sadj scaled vector)
Nunique = length(Sadj);
% NFFT is length of full spectrum (2 * size of spectrum w/o DC)
NFFT = 2*(Nunique - 1);

% build complex spectrum
Scomplex = complex(Sadj.*cos(Sphase), Sadj.*sin(Sphase));

% assign indices into Sfull for the two "sections"
indx{1} = 1:Nunique;
% second portion
indx{2} = (Nunique+1):NFFT;

% assign to Sfull
Sfull = zeros(1, NFFT);
Sfull(indx{1}) = Scomplex;
% second section is computed as:
%	(1) take fftred(1:(end-1)), since final point (fftred(end)) 
% 		 is common to both sections
% 	(2) flip the fftred section around using fliplr (reverse order)
% 	(3) take complex conjugate of flipped fftred
Sfull(indx{2}) = conj(fliplr(Scomplex(2:(end-1))));

% take ifft
ssyn = ifft(Sfull);

subplot(313)
plot(abs(Sfull), '.-')
grid
title('abs(Sfull)')


figure(2)

subplot(311)
plot(s)
grid
title('s original')

subplot(312)
plot(real(ssyn(1:N)))
grid
title('s synthesized')

serr = s - real(ssyn(1:N));
fprintf('max serr = %e\n\n', max(abs(serr)));
subplot(313)
plot(serr)
grid
title('s - ssn')

%------------------------------------------------------------------------------
%------------------------------------------------------------------------------
%------------------------------------------------------------------------------
% test synverse function
%------------------------------------------------------------------------------
%------------------------------------------------------------------------------
[sv, SVfull, sim]  = synverse(Smag, Sphase, 'DC', 'yes', 'Scale', N);

figure(3)
plot(sv)
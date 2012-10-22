function [sadj, Sfull, Magnorm, f] = compensate_signal(s, calfreq, calmag, Fs, corr_frange, varargin)
%------------------------------------------------------------------------
% [sadj, Sfull, Magnorm, f] = compensate_signal(s, calfreq, calmag, Fs, corr_frange)
%------------------------------------------------------------------------
% 
% Description
% 
%------------------------------------------------------------------------
% Input Arguments:
% 	
% 
% 	Options:
% 		Method
% 		
% 		Normalize
% 		
% 		Lowcut
% 		
% 		Level (only applicable for COMPRESS method!)
% 		
% 
% Output Arguments:
% 	Output	output info
%
%------------------------------------------------------------------------
% See also: NICal
%------------------------------------------------------------------------

%------------------------------------------------------------------------
% Sharad J. Shanbhag
% sshanbhag@neomed.edu
%------------------------------------------------------------------------
% Created: XX XXXX, 2011 (SJS)
%
% Revisions:
%	1 Oct 2012 (SJS): working on method # 2
%	8 Oct 2012 (SJS): implemented COMPRESS method
%	9 Oct 2012 (SJS): added LEVEL option to specify target level
%------------------------------------------------------------------------
% TO DO:
%------------------------------------------------------------------------

%------------------------------------------------------------------------
%------------------------------------------------------------------------
% define some constants
%------------------------------------------------------------------------
%------------------------------------------------------------------------
% arbitrary minimum dB value
MIN_DB = -120;
% need to have a small, but non-zero value when taking log, so set that here
ZERO_VAL = 1e-17;

%------------------------------------------------------------------------
%------------------------------------------------------------------------
% define defaults for settings
%------------------------------------------------------------------------
%------------------------------------------------------------------------
% default method
COMPMETHOD = 'BOOST';
% default low frequency cutoff.  4 kHz is used to protect speakers!
LOWCUT = 4000;
% normalize
NORMALIZE = 0;
% level sets target flattening level; if 0 (default), will set from cal data 
% depending on method
LEVEL = 0;

%------------------------------------------------------------------------
%------------------------------------------------------------------------
% check input arguments
%------------------------------------------------------------------------
%------------------------------------------------------------------------
nvararg = length(varargin);
if nvararg
	aindex = 1;
	while aindex <= nvararg
		switch(upper(varargin{aindex}))
			
			% select method
			case 'METHOD'
				mtype = upper(varargin{aindex + 1});
				switch mtype
					case 'BOOST'
						COMPMETHOD = 'BOOST';
					case 'ATTEN'
						COMPMETHOD = 'ATTEN';
					case 'COMPRESS'
						COMPMETHOD = 'COMPRESS';
					otherwise
						fprintf('%s: unknown compensation method %s\n', mfilename, mtype);
						fprintf('\tUsing default, BOOST method\n');
						COMPMETHOD = 'BOOST';
				end
				aindex = aindex + 2;
				clear mtype;
			
			% set LOWCUT (low frequency cutoff) option & frequency
			case 'LOWCUT'
				lval = varargin{aindex + 1};
				if strcmpi(lval, 'OFF')
					LOWCUT = 0;
				elseif isnumeric(lval)
					LOWCUT = lval;
				end
				aindex = aindex + 2;
				clear lval;
				
			% set normalization of output signal
			case 'NORMALIZE'
				nval = varargin{aindex + 1};
				if strcmpi(nval, 'ON')
					NORMALIZE = 1;
				elseif isnumeric(nval)
					NORMALIZE = nval;
				else
					NORMALIZE = 0;
				end
				aindex = aindex + 2;
				clear nval;
				
			% set level
			case 'LEVEL'
				lval = varargin{aindex + 1};
				if isnumeric(lval)
					if lval <= 0
						error('%s: LEVEL value must be greater than zero!', mfilename);
					else
						LEVEL = lval;
					end
				else
					error('%s: invalid LEVEL value (%s)', mfilename, lval);
				end
				aindex = aindex + 2;
				clear lval;
				
			otherwise
				error('%s: Unknown option %s', mfilename, varargin{aindex});
		end		% END SWITCH
	end		% END WHILE aindex
end		% END IF nvararg


%------------------------------------------------------------------------
%------------------------------------------------------------------------
% get spectrum of raw signal
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
% apply correction using BOOST method
%------------------------------------------------------------------------
%------------------------------------------------------------------------
% procedure:	find additive compensation values for frequency range 
%					for which there are calibration data and apply to FFT, 
%					then iFFT to get corrected version
%------------------------------------------------------------------------
%------------------------------------------------------------------------
if strcmpi(COMPMETHOD, 'BOOST')
	% normalize to peak of xfer function
	% find peak magnitude
	peakmag = max(calmag);
	% normalize by finding deviation from peak
	Magnorm = peakmag - calmag;

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
end


%------------------------------------------------------------------------
%------------------------------------------------------------------------
% apply correction using ATTEN method
%------------------------------------------------------------------------
%------------------------------------------------------------------------
% procedure:	find subtractive compensation values for frequency range 
%					for which there are calibration data and apply to FFT, 
%					then iFFT to get corrected version
%------------------------------------------------------------------------
%------------------------------------------------------------------------
if strcmpi(COMPMETHOD, 'ATTEN')
	% normalize to peak of xfer function
	% find peak magnitude
	minmag = min(calmag);
	% normalize by finding deviation from peak
	Magnorm = minmag - calmag;
	
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
end

%------------------------------------------------------------------------
%------------------------------------------------------------------------
% apply correction using COMPRESS method
%------------------------------------------------------------------------
%------------------------------------------------------------------------
% procedure:	find subtractive compensation values for frequency range 
%					for which there are calibration data and apply to FFT, 
%					then iFFT to get corrected version
%------------------------------------------------------------------------
%------------------------------------------------------------------------
% some assumptions:
% 					magnitude values are in ACTUAL, dB SPL range.  
% 						¡this algorithm blows up for negative magnitudes!
%------------------------------------------------------------------------
%------------------------------------------------------------------------
if strcmpi(COMPMETHOD, 'COMPRESS')
	
	% check if LEVEL was specified
	if LEVEL
		% normalize by finding deviation from peak
		Magnorm = LEVEL - calmag;
	else
		% find max and min in magnitude spectrum
		maxmag = max(calmag);
		minmag = min(calmag);
		% compute middle value
		midmag = ((maxmag - minmag) / 2) + minmag;
		% normalize by finding deviation from peak
		Magnorm = midmag - calmag;
	end
	
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
end

%------------------------------------------------------------------------
%------------------------------------------------------------------------
% normalize output if desired
%------------------------------------------------------------------------
%------------------------------------------------------------------------
if NORMALIZE
	keyboard
	sadj = NORMALIZE * normalize(sadj);
end


function [y, start_bins, end_bins] = block_rms(x, rmsbins)
%--------------------------------------------------------------------------
% [y, start_bins, end_bins] = block_rms(x, rmsbins)
%--------------------------------------------------------------------------
%
% computes root-mean-square (rms) of vector x in windows of rmsbins size
%
%--------------------------------------------------------------------------
% Input Arguments:
%	x			signal (vector)
%	rmsbins	size (in samples) of rms window
%
% Output Arguments:
% 	y			rms values
% 	start_bins	rms window start locations
% 	end_bins		rms window end locations
%--------------------------------------------------------------------------
% See Also: rms, dbspl
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% Sharad J. Shanbhag
% sshanbhag@neomed.edu
%--------------------------------------------------------------------------
% Created:	????
%
% Revision History:
%	4 Apr 2017 (SJS): added comments, moved copy to GeneralUtilities 
%		toolbox
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

% compute starting and end locations of rms windows
start_bins = 1:rmsbins:length(x);
end_bins = start_bins + (rmsbins - 1);
% keep only valid indices
validindx = end_bins <= length(x);
start_bins = start_bins(validindx);
end_bins = end_bins(validindx);

% perform rms on chunks of y vector
nrms = length(start_bins);
y = zeros(nrms, 1);
for n = 1:nrms
	y(n) = rms(x(start_bins(n):end_bins(n)));
end
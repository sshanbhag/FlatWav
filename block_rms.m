function [y, start_bins, end_bins] = block_rms(x, rmsbins)
% computes rms of vector x in windows of rmsbins size

start_bins = 1:rmsbins:length(x);
end_bins = start_bins + (rmsbins - 1);
validindx = end_bins <= length(x);

start_bins = start_bins(validindx);
end_bins = end_bins(validindx);

nrms = length(start_bins);
y = zeros(nrms, 1);
for n = 1:nrms
	y(n) = rms(x(start_bins(n):end_bins(n)));
end
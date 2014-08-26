% test_rmsplot.m
load rmsvals.mat
% test rms sig plot

if any(	[	isempty(respFs) ...
				isempty(rawresp)])
	error('No data for db analysis');
end
	
PeakRMSWindow = 5;
rmsbins = ms2samples(PeakRMSWindow, respFs)

fprintf('Creating dBFigure\n\n');
dBFigure = figure(10);
% 	set(handles.dBFigure, 'Position', [10 235 973 500]);
% create subplot axes
P.rawdb = subplot(211);
P.adjdb = subplot(212);

[rawrms, start_bins, end_bins] = block_rms(rawresp, rmsbins);

% start_bins = 1:rmsbins:length(rawresp);
% end_bins = start_bins + (rmsbins - 1);
% validindx = end_bins <= length(rawresp);
% 
% start_bins = start_bins(validindx);
% end_bins = end_bins(validindx);
% 
% nrms = length(start_bins);
% rmsvals = zeros(nrms, 1);
% for n = 1:nrms
% 	rmsvals(n) = rms(rawresp(start_bins(n):end_bins(n)));
% end

nrms = length(rawrms)
x = zeros(2*nrms, 1);
y = zeros(2*nrms, 1);
for n = 1:nrms
	x(2 * (n - 1) + 1) = start_bins(n);
	x(2 * (n - 1) + 2) = end_bins(n);
	y(2 * (n - 1) + 1) = rawrms(n);
	y(2 * (n - 1) + 2) = rawrms(n);
end
y = db(y);
plotyy(bin2ms(1:length(rawresp), respFs), rawresp, ...
			bin2ms(x, respFs), y + max(abs(y)))

return
%




% find peak and peak index of rms values
[rawrespmax.val, rawrespmax.indx] = max(rawrms);
% compute peak dB SPL
rawdBSPL = dbspl(rawrespmax.val);

% find max point (in milliseconds)
xval = rmsbins * rawrespmax.indx + (rmsbins ./ 2);
xval = fix(bin2ms(xval, respFs));

% plot
tvec = 1000 * (0:(length(rawresp)-1)) ./ respFs;
subplot(P.rawdb)
x = bin2ms(rmsbins * (1:(length(rawrms))), respFs);
plotyy(P.rawdb, tvec, rawresp,  ...
			x, dbspl(rawrms));

title(P.rawdb, 'Signal (V)')
ylabel(P.rawdb, 'Raw', 'Color', 'b')
set(P.rawdb, 'XTickLabel', []);
xlim(P.rawdb, [min(tvec) max(tvec)])
% get ticks
time_ticks = get(P.rawdb, 'XTick');	



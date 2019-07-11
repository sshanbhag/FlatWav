function [caldata, varargout] = padCalData
%------------------------------------------------------------------------
% 
%------------------------------------------------------------------------
% FlatCal
%------------------------------------------------------------------------
% kludge to pad calibration data beyond calibration range
%------------------------------------------------------------------------
% Input Arguments:
% 
% Output Arguments:
%
%------------------------------------------------------------------------
% See also: FlatWav
%------------------------------------------------------------------------

%------------------------------------------------------------------------
% Sharad J. Shanbhag
% sshanbhag@neomed.edu
%------------------------------------------------------------------------
% Created: 11 July, 2019 (SJS)
%
%------------------------------------------------------------------------
% TO DO:
%------------------------------------------------------------------------

%% load caldata
[calfile, calpath] = uigetfile( {'*.cal'; '*_cal.mat'}, ...
									'Load calibration data from file...');
if calfile ~=0
	datafile = fullfile(calpath, calfile);
	oCal = load_headphone_cal(datafile);
else
	fprintf('Cancelled\n');
	return
end

%%
% common format fields
array_fields = {'mag', ...
'phase', ...
'dist', ...
'mag_stderr', ...
'phase_stderr', ...
'background', ...
'background_stderr', ...
'dist_stderr', ...
'phase_us', ...
'maginv'};
cell_fields = {'atten', ...
					'magsraw', ...
					'magsV'};
%%
figure;
oH = gca;
plot(	oH, ...
		0.001*oCal.freq, ...
		oCal.mag(1, :), 'o-');
ylim(oH, ...
		[0.9*min(oCal.mag(1, :)) 1.1*max(oCal.mag(1, :))]);
grid(	oH, 'on');
grid minor
xlabel('kHz');
ylabel('dB SPL');
title(calfile, 'Interpreter', 'none');
%% figure out freqs
if oCal.settings.Fmin <= oCal.settings.Fstep
	error('%s: calibration Fmin already lower than Fstep', mfilename);
end

% determine pad freqs
fpad = 0:oCal.settings.Fstep:(oCal.settings.Fmin - oCal.settings.Fstep);
padlen = length(fpad);
if padlen == 0
	error('not enough values to pad');
end

%% make copy and assign padded values
nCal = oCal;

% freq is easy to do - 
nCal.freq = [fpad oCal.freq];
nCal.range(1) = fpad(1);

% loop through cell arrays
for N = 1:length(cell_fields)
	for c = 1:length(oCal.(cell_fields{N}))
		% get size of array
		[nr, nc] = size(oCal.(cell_fields{N})); %#ok<ASGLU>
		% get pad value
		pval = oCal.(cell_fields{N}){c}(1, :);
		% create pad array
		parr = pval .* ones(padlen, nc);
		% append pad to array, and assign to new cal
		nCal.(cell_fields{N}){c} = [parr; oCal.(cell_fields{N}){c}];
	end
end

%%
% loop through regular arrays
for N = 1:length(array_fields)
	% get pad value
	pval = oCal.(array_fields{N})(:, 1);
	% create pad array
	parr = pval .* ones(2, padlen);
	% append pad to array, assign to new cal
	nCal.(array_fields{N}) = [parr oCal.(array_fields{N})];
end

%% plot padded data
hold on
plot(	oH, ...
		0.001*nCal.freq, ...
		nCal.mag(1, :), 'r.:');
legend({'original', 'padded'});
hold off

%%  save padded cal file
% get base of calfile
[~, basestr, extstr] = fileparts(calfile);
% append _pad to name
tmpfile = [basestr '_PAD' extstr];

[newfile, newpath] = uiputfile( fullfile(calpath, tmpfile), ...
									'Write padded calibration data to file...');
if newfile ~=0	
	% remove unnecessary fields
	caldata = rmfield(nCal, {'phase_us', 'maginv'});
	save(fullfile(newpath, newfile), 'caldata', '-MAT');
else
	fprintf('Cancelled\n');
	return
end

%% output
if nargout > 1
	varargout{1} = load_headphone_cal(fullfile(newpath, newfile));
end
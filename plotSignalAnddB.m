function varargout = plotSignalAnddB(signal, rmswin, Fs, varargin)
%-------------------------------------------------------------------------
% [dBFigure, dBAxes, out] = plotSignalAnddB(signal, rmswin, ...
% 																			Fs, 'axes', axH)
%-------------------------------------------------------------------------
%  TytoLogy:PlotTools
%-------------------------------------------------------------------------
% Plots signal along with dbSPL in windows
%------------------------------------------------------------------------
% Input Arguments:
%	signal		vector of data to plot
% 	rmswin		time window (milliseconds) for computing RMS and dB SPL
% 	Fs				sampling rate for signal (samples/s)
% 	
% 	Optional:
% 		'axes'			<axis handle>	specify axis handle for plotting
% 		'dBSPL'			<Volts to Pa conversion factor>
% 		'SIGNALNAME'	<signal trace name>
% 		'SIGNALCOLOR'			<signal trace color>
% 		'SIGNALSTYLE'			<signal trace style>
% 		'DBCOLOR'				<dB trace color>
% 		'DBSTYLE'				<db trace style>
% 		'DBLINEWIDTH'			<db trace line width>
% 		'DBMARKER'				<db peak marker symbol>
% 		'DBMARKERSIZE'			<db peak marker size (points)>
% 		'DBMARKERCOLOR'		<db peak marker size (points)>
%
% Output Arguments:
%	dBFigure		figure handle of plot
%	dbAxex		axis handle of plot
%	out			output struct with fields:
% 		dB_max		peak dB (or dB SPL value if dBSPL option set)
% 		rms_max		peak RMS value = maxval;
% 		max_bin		index of RMS window for peak dB
% 		max_time		time of peak dB window
% 		rms_vals		vector of rms values
% 		dBSPL			vector of dB SPL values (if dBSPL option set)
% 		dB				vector of dB values (if dBSPL conversion not provided)
%	   
%------------------------------------------------------------------------
% See Also: plot, dbspl, db
%------------------------------------------------------------------------

%-------------------------------------------------------------------------
% Sharad J Shanbhag
% sshanbhag@neomed.edu
%-------------------------------------------------------------------------
% Created (as separate function): 30 September 2016 (SJS)
% 				pulled out of FlatWav.m
% 
% Revisions:
%	6 Jun 2019 (SJS):
%		- updated comments
%		- added out struct
%		- moving to PlotTools
%-------------------------------------------------------------------------

%-------------------------------------------------------------------------
% definitions/defaults
%-------------------------------------------------------------------------
VtoPa = 0;
rmswin = 5;
sigName = '';
sigColor = 'b';
sigStyle = '-';
dBColor = 'k';
dBStyle = '-';
dBLineWidth = 2;
dBMarker = '*';
dBMarkerSize = 10;
dBMarkerColor = dBColor;

%-------------------------------------------------------------------------
% check input arguments
%-------------------------------------------------------------------------
nvararg = length(varargin);
if nvararg
	aindex = 1;
	while aindex <= nvararg
		switch(upper(varargin{aindex}))
			%-----------------------------
			% select axes
			%-----------------------------
			case 'AXES'
				if ~ishandle(varargin{aindex + 1})
					error('%s: invalid axes handle', mfilename);
				end
				dBAxes = varargin{aindex + 1};
				aindex = aindex + 2;
			%-----------------------------
			% set dBSPL option
			%-----------------------------
			case 'DBSPL'
				% store conversion factor
				VtoPa = varargin{aindex + 1};
				aindex = aindex + 2;
			%-----------------------------
			% set signal trace name
			%-----------------------------
			case 'SIGNALNAME'
				sigName = varargin{aindex + 1};
				aindex = aindex + 2;
			%-----------------------------
			% set signal trace color
			%-----------------------------
			case 'SIGNALCOLOR'
				sigColor = varargin{aindex + 1};
				aindex = aindex + 2;				
			%-----------------------------
			% set signal trace style
			%-----------------------------
			case 'SIGNALSTYLE'
				sigStyle = varargin{aindex + 1};
				aindex = aindex + 2;				
			%-----------------------------
			% set dB trace color
			%-----------------------------
			case 'DBCOLOR'
				dbColor = varargin{aindex + 1}; %#ok<NASGU>
				aindex = aindex + 2;				
			%-----------------------------
			% set db trace style
			%-----------------------------
			case 'DBSTYLE'
				dbStyle = varargin{aindex + 1}; %#ok<NASGU>
				aindex = aindex + 2;				
			%-----------------------------
			% set db trace line width
			%-----------------------------
			case 'DBLINEWIDTH'
				dBLineWidth = varargin{aindex + 1};
				aindex = aindex + 2;				
			%-----------------------------
			% set db peak marker symbol
			%-----------------------------
			case 'DBMARKER'
				dBMarker = varargin{aindex + 1};
				aindex = aindex + 2;				
			%-----------------------------
			% set db peak marker size (points)
			%-----------------------------
			case 'DBMARKERSIZE'
				dBMarkerSize = varargin{aindex + 1};
				aindex = aindex + 2;
			%-----------------------------
			% set db peak marker size (points)
			%-----------------------------
			case 'DBMARKERCOLOR'
				dBMarkerColor = varargin{aindex + 1};
				aindex = aindex + 2;				
			%-----------------------------
			% invalid argument
			%-----------------------------
			otherwise
				error('%s: Unknown option %s', mfilename, varargin{aindex});
		end		% END SWITCH
	end		% END WHILE aindex
else
	dBFigure = figure;
	dBAxes = axes;	
end		% END IF nvararg

%-------------------------------------------------------------------------
% process inputs
%-------------------------------------------------------------------------
% convert rmswindow from milliseconds to # of bins
rmsbins = ms2samples(rmswin, Fs);
% compute rms of  response, plot
[rawrms, startbins, endbins] = block_rms(signal, rmsbins);
% find peak and peak index of rms values
[maxval, maxindx] = max(rawrms);

if VtoPa
	% compute peak dB SPL
	sigdBSPL = dbspl(VtoPa*maxval);
	% display value
	dbtext = sprintf('%s Peak dB SPL: %.2f\n', sigName, sigdBSPL);
	fprintf('%s\n', dbtext);
else
	% just use dB
	sigdBSPL = db(maxval);
	% display value
	dbtext = sprintf('%s Peak dB: %.2f\n', sigName, sigdBSPL);
	fprintf('%s\n', dbtext);
end

% find max point (in milliseconds)
xval = rmsbins * maxindx - (rmsbins ./ 2);
xval = fix(bin2ms(xval, Fs));

% build trace for dB data
nrms = length(rawrms);
x = zeros(2*nrms, 1);
y = zeros(2*nrms, 1);
for n = 1:nrms
	x(2 * (n - 1) + 1) = startbins(n);
	x(2 * (n - 1) + 2) = endbins(n);
	y(2 * (n - 1) + 1) = rawrms(n);
	y(2 * (n - 1) + 2) = rawrms(n);
end
x = bin2ms(x, Fs);
if VtoPa
	y = dbspl(VtoPa*y);
else
	y = db(y);
end

% response data
tvec = bin2ms( (1:length(signal))-1, Fs);
yresp = max(y) * normalize(signal);
xlimits = [min(tvec) max(tvec)];
ylimits = [min(yresp) 1.05*max(y)];

%-------------------------------------------------------------------------
% Plot!
%-------------------------------------------------------------------------

% plot signal data
plot(dBAxes, tvec, yresp, [sigColor sigStyle]);

hold(dBAxes, 'on')
	% plot dB data
	plot(dBAxes, x, y, [dBColor dBStyle], 'LineWidth', dBLineWidth);
	% plot dB peak
	plot(dBAxes, xval, ylimits(2), [dBMarkerColor dBMarker], ...
												'MarkerSize', dBMarkerSize);
hold(dBAxes, 'off')
ylabel(dBAxes, sigName, 'Color', sigColor)
xlim(dBAxes, xlimits);
ylim(dBAxes, ylimits);
grid(dBAxes, 'on');
th = text(xval, 1.05*ylimits(2), sprintf('  %.2f', sigdBSPL), ...
																'Parent', dBAxes);
set(th,	'FontSize', 12, ...
			'FontWeight', 'bold', ...
			'Color', sigColor, ...
			'Interpreter', 'none');

%-------------------------------------------------------------------------
% assign outputs
%-------------------------------------------------------------------------
if nargout > 0
	varargout{1} = dBFigure;
	varargout{2} = dBAxes;

	out.dB_max = sigdBSPL;
	out.rms_max = maxval;
	out.max_bin = maxindx;
	out.max_time = xval;
	out.rms_vals = raw_rms;
	if VtoPa
		% dbSPL vals
		out.dBSPL = dbspl(VtoPa*raw_rms);
		out.dB = [];
	else
		% just use dB
		out.dB = db(VtoPa*raw_rms);
		out.dBSPL = [];
	end
	varargout{3} = out;
end

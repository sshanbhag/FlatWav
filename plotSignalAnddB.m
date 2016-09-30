function varargout = plotSignalAnddB(signal, rmswin, Fs, varargin)
%-------------------------------------------------------------------------
% varargout = plotSignalAnddB(signal, rmswin, Fs, varargin)
%-------------------------------------------------------------------------
% Plots signal along with dbSPL in windows
%------------------------------------------------------------------------
% Input Arguments:
%	signal
% 	rmswin
% 	Fs
% 	
% 	Optional:
% 		'axes'
% 		
% Output Arguments:
%   
%------------------------------------------------------------------------

%-------------------------------------------------------------------------
% Sharad J Shanbhag
% sshanbhag@neomed.edu
%-------------------------------------------------------------------------
% Created (as separate function): 30 September 2016 (SJS)
% 				pulled out of FlatWav.m
% 
% Revisions:
%-------------------------------------------------------------------------

%-------------------------------------------------------------------------
% definitions
%-------------------------------------------------------------------------
VtoPa = 0;
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
if nargout
	varargout{1} = dBFigure;
	varargout{2} = dBAxes;
end

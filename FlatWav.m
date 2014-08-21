function varargout = FlatWav(varargin)
% FLATWAV MATLAB code for FlatWav.fig
%      FLATWAV, by itself, creates a new FLATWAV or raises the existing
%      singleton*.
%
%      H = FLATWAV returns the handle to a new FLATWAV or the handle to
%      the existing singleton*.
%
%      FLATWAV('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in FLATWAV.M with the given input arguments.
%
%      FLATWAV('Property','Value',...) creates a new FLATWAV or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before FlatWav_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to FlatWav_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help FlatWav

% Last Modified by GUIDE v2.5 21-Aug-2014 16:40:14

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @FlatWav_OpeningFcn, ...
                   'gui_OutputFcn',  @FlatWav_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT

%******************************************************************************
%******************************************************************************
%******************************************************************************
% Essential Functions
%******************************************************************************
%******************************************************************************
%******************************************************************************

%------------------------------------------------------------------------------
%------------------------------------------------------------------------------
% --- Executes just before FlatWav is made visible.
%------------------------------------------------------------------------------
function FlatWav_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to FlatWav (see VARARGIN)

	% Choose default command line output for FlatWav
	handles.output = hObject;
	% Update handles structure
	guidata(hObject, handles);

	%----------------------------------------------------------
	%----------------------------------------------------------
	% Setup Paths
	%----------------------------------------------------------
	%----------------------------------------------------------
	disp([mfilename ': checking paths'])
	if ispc
		% directory when using installed version:
		%	pdir = ['C:\TytoLogy\TytoLogySettings\' getenv('USERNAME')];
		% development tree
		pdir = ['C:\Users\sshanbhag\Code\Matlab\TytoLogy\TytoLogySettings\' getenv('USERNAME')];
	else ismac
		pdir = ['~/Work/Code/Matlab/dev/TytoLogy/TytoLogySettings/' getenv('USER')];
	end
			
	if isempty(which('ms2samples'))
		% could not find the RPload.m function (which is in TytoLogy
		% toolbox) which suggests that the paths are not set or are 
		% incorrect for this setup.  load the paths using the tytopaths program.
		run(fullfile(pdir, 'tytopaths'))
	else
		disp([mfilename ': paths ok, launching programn'])
	end

	%--------------------------------------------------
	%--------------------------------------------------
	% SYNTH SETTINGS
	%--------------------------------------------------
	%--------------------------------------------------
	% set output signal to synth (vs. wav) and update GUI
	handles.SignalMode = 'SYNTH';
	update_ui_val(handles.SynthSignalButton, 1);
	update_ui_val(handles.WavSignalButton, 0);
	% initialize S synth parameter structure
	handles.S = FlatWav_buildS;
	guidata(hObject, handles);
	% initialize tone, noise, sweep structs
	typenum = 1;
	handles.Tone.type = handles.S.Types{typenum};
	for n = 1:handles.S.Nparam(typenum);
		handles.Tone.(handles.S.Param{typenum}{n}) = handles.S.DefaultVals{typenum}(n);
	end
	typenum = 2;
	handles.Noise.type = handles.S.Types{typenum};
	for n = 1:handles.S.Nparam(typenum);
		handles.Noise.(handles.S.Param{typenum}{n}) = handles.S.DefaultVals{typenum}(n);
	end
	typenum = 3;
	handles.Sweep.type = handles.S.Types{typenum};
	for n = 1:handles.S.Nparam(typenum);
		handles.Sweep.(handles.S.Param{typenum}{n}) = handles.S.DefaultVals{typenum}(n);
	end
	% set current synth object to Noise
	handles.synth = handles.Noise;
	handles.SynthIndex = 2;
	handles.SynthType = handles.S.Types{handles.SynthIndex};
	% set Analysis window
	handles.Awindow = [0 handles.S.DefaultVals{typenum}(1)];

	%--------------------------------------------------
	%--------------------------------------------------
	% update GUI and synth
	%--------------------------------------------------
	%--------------------------------------------------
	guidata(hObject, handles);
	updateGuiFromSynth(hObject, handles)
	guidata(hObject, handles);
	updateSynthFromGui(hObject, handles);
	guidata(hObject, handles);
	% set Raw and Adj dB text to default vals and hide them
	update_ui_str(handles.RawdBText, sprintf('Raw dB SPL: ---'));
	update_ui_str(handles.AdjdBText, sprintf('Adj dB SPL: ---'));
	hide_uictrl(handles.RawdBText);
	hide_uictrl(handles.AdjdBText);
	% assign blank handle to PlaySignalFig
	handles.PlaySignalFig = [];
	% signal IO figure
	handles.IOfigure = [];
	guidata(hObject, handles);
	
	
	%--------------------------------------------------
	%--------------------------------------------------
	% COMPENSATION SETTINGS
	%--------------------------------------------------
	%--------------------------------------------------
	% reset string in CompMethodCtrl
	set(handles.CompMethodCtrl, 'string', 'none|normalize|atten|boost|compress');
	% set compensation method to 1 ('atten')
	handles.CompMethod = 1;
	update_ui_val(handles.CompMethodCtrl, handles.CompMethod);
	% set correction freq range (in Hz) and update ctrl values
	handles.CorrFrange = [100 10000];
	update_ui_str(handles.CorrFminCtrl, handles.CorrFrange(1));
	update_ui_str(handles.CorrFmaxCtrl, handles.CorrFrange(2));
	% default normalize status
	handles.Normalize = 'on';
	handles.NormalizeValue = 1.0;
	update_ui_val(handles.NormalizeCtrl, 1);
	update_ui_str(handles.NormalizePeakCtrl, handles.NormalizeValue);
	show_uictrl(handles.NormalizePeakCtrl);
	show_uictrl(handles.NormalizePeakText);
	% default LowCut options
	handles.LowCut = 'off';
	handles.LowCutFreq = read_ui_str(handles.LowCutFreqCtrl, 'n');
	if strcmpi(handles.LowCut, 'off')
		disable_ui(handles.LowCutFreqText);
		disable_ui(handles.LowCutFreqCtrl);
	end
	% target SPL
	handles.TargetSPL = 65;
	update_ui_str(handles.TargetSPLCtrl, handles.TargetSPL);
	guidata(hObject, handles);

	%--------------------------------------------------
	%--------------------------------------------------
	% spectrum settings
	%--------------------------------------------------
	%--------------------------------------------------
	handles.SpectrumWindow = 512;
	handles.ColorMap = 'gray';
	guidata(hObject, handles);
	
	%--------------------------------------------------
	%--------------------------------------------------
	% set initial state for sounds
	%--------------------------------------------------
	% wavdata struct holds information about wav file.
	% create blank wavdata struct, update gui
	handles.wavdata = struct(	'datafile', '', 'raw', [], 'fs', [], ...
										'nbits', [], 'opts', []);
	update_ui_str(handles.WavFilenameCtrl, '');
	update_ui_str(handles.WaveInfoCtrl, 'no wav loaded');
	% create empty raw and adj vectors
	handles.raw = [];
	handles.adj = [];
	guidata(hObject, handles);

	%--------------------------------------------------
	%--------------------------------------------------
	% fake cal data
	%--------------------------------------------------
	handles.cal = fake_caldata('freqs', 1:10:(handles.S.Fs / 2));
	handles.cal.mag = 90 * handles.cal.mag;
	guidata(hObject, handles);
	plot(handles.CalibrationAxes, 0.001*handles.cal.freq, handles.cal.mag(1, :), '.-');
	ylim([0 100]);
	
	%--------------------------------------------------
	%--------------------------------------------------
	% calibration processing settings
	%--------------------------------------------------
	handles.SmoothMethod = 1;
	handles.SmoothVal1 = 3;
	handles.SmoothVal2 = 5;
	guidata(hObject, handles)
	SmoothCalCtrl_Callback(hObject, eventdata, handles);

	%--------------------------------------------------
	%--------------------------------------------------
	% output settings
	%--------------------------------------------------
	% options for Output Device are 'winsound' and 'NI-DAQ'
	handles.OutputDevice = 'winsound';
	handles.IOpause = 0.5;
	guidata(hObject, handles);
	
	%--------------------------------------------------
	%--------------------------------------------------
	% Mic settings
	%--------------------------------------------------
	handles.MicSensitivity = 1;
	handles.MicGaindB = 0;
	handles.MicGain = invdb(handles.MicGaindB);
	handles.VtoPa = (1/handles.MicGain) * (1/handles.MicSensitivity);
	guidata(hObject, handles);

	
%******************************************************************************
%******************************************************************************
%******************************************************************************


%******************************************************************************
%******************************************************************************
%******************************************************************************
%------------------------------------------------------------------------------
% --- Outputs from this function are returned to the command line.
%------------------------------------------------------------------------------
function varargout = FlatWav_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

	% Get default command line output from handles structure
	varargout{1} = handles.output;
%------------------------------------------------------------------------------
%******************************************************************************
%******************************************************************************
%******************************************************************************


%******************************************************************************
%******************************************************************************
%******************************************************************************
% ACTION BUTTON CONTROL Callbacks
%******************************************************************************
%******************************************************************************
%******************************************************************************

%------------------------------------------------------------------------------
% --- Executes on button press in UpdateSignalCtrl.
%------------------------------------------------------------------------------
function UpdateSignalCtrl_Callback(hObject, eventdata, handles)
%------------------------------------------------------------------------------
% using settings from GUI, update the test signal, compensated signal, and
% plots of the signals
%------------------------------------------------------------------------------
	disp('updating signal display...')
	%--------------------------------------------------------------------
	% Signal Mode is either WAV (.wav file) or SYNTH (synthesized)
	%--------------------------------------------------------------------
	switch handles.SignalMode
		
		case 'WAV'
			% load wav file
			[wavdata.raw, wavdata.Fs, wavdata.nbits, wavdata.opts] = ...
														wavread(handles.wavdata.datafile);
			% apply ramp (short, just to ensure zeros at beginning and end of
			% stimulus)
			wavdata.datafile = handles.wavdata.datafile;
			wavdata.raw = sin2array(wavdata.raw', 0.5, wavdata.Fs);
			% store raw vector in handles.
			handles.raw = wavdata.raw;
			handles.S.Fs = wavdata.Fs;
			handles.wavdata = wavdata;
			update_ui_str(handles.FsCtrl, handles.S.Fs);
			
		case 'SYNTH'
			% update Synth data (handles.synth) from GUI
			updateSynthFromGui(hObject, handles);
			% store changes in handles
			guidata(hObject, handles);
			% local copy of synth for sake of brevity
			synth = handles.synth;
			% act according to type of synth signal
			switch synth.type
				case 'tone'
					% create tone
					handles.raw = synmonosine(	synth.dur, ...
														handles.S.Fs, ...
														synth.freq, ...
														synth.amp, 0);

				case 'noise'
					% create noise
					handles.raw = synmononoise_fft(	synth.dur, ...
																handles.S.Fs, ...
																synth.fmin, ...
																synth.fmax, ...
																synth.amp, 0);
					% kludge to scale amplitude properly
					handles.raw = synth.amp * normalize(handles.raw);
					
				case 'sweep'
					% create FM sweep (via wrapper around chirp() matlab
					% function)
					handles.raw = synmonosweep(	synth.dur, ...
															handles.S.Fs, ...
															synth.fmin, ...
															synth.fmax, ...
															synth.amp, 0);
			end	% END switch synth.type
			
	end	% END switch handles.SignalMode
	
	% update analysis window
	% check if analysis window is beyond length of signal
	if  ms2samples(handles.Awindow(2), handles.S.Fs) > length(handles.raw)
		% if so, reset to duration of signal
		handles.Awindow(2) = floor(bin2ms(length(handles.raw), handles.S.Fs));
		update_ui_str(handles.AnalysisEndCtrl, handles.Awindow(2));
		guidata(hObject, handles);
	end
	% find bins for analysis
	bin = ms2samples(handles.Awindow, handles.S.Fs);
	if bin(1) == 0
		bin(1) = 1;
	end
	
	% take fft of raw data
	[handles.fraw, handles.magraw, handles.phiraw] = ...
									daqdbfullfft(	handles.raw(bin(1):bin(2)), ...
														handles.S.Fs, ...
														length(handles.raw(bin(1):bin(2))));
	guidata(hObject, handles);

	%--------------------------------------------------------------------
	% apply compensation method
	% CompMethod is value of 1, 2, 3, 4 (value of CompMethodCtrl) which 
	% is a pull-down menu.  Map these to string values
	%--------------------------------------------------------------------
	switch handles.CompMethod
		case 1
			method = 'NONE';
		case 2
			method = 'NORMALIZE';
		case 3
			method = 'ATTEN';
		case 4
			method = 'BOOST';
		case 5
			method = 'COMPRESS';
	end
	
	% check low freq cutoff setting
	if strcmp(handles.LowCut, 'off')
		lowcut = 'off';
	else
		lowcut = handles.LowCutFreq;
	end
	
	if strcmpi(method, 'NONE')
		handles.adj = handles.raw;
	
	elseif strcmpi(method, 'NORMALIZE')
		% normalize only
		handles.adj = handles.NormalizeValue * normalize(handles.raw);
		
	elseif strcmpi(handles.Normalize, 'off')
		handles.adj = compensate_signal(	handles.raw, ...
											handles.cal.freq, ...
											handles.cal.mag(1, :), ...
											handles.S.Fs, ...
											handles.CorrFrange, ...
											'Method', method, ...
											'Normalize', 'off', ...
											'Lowcut', lowcut, ...
											'Level', handles.TargetSPL);
	else
		handles.adj = compensate_signal(	handles.raw, ...
											handles.cal.freq, ...
											handles.cal.mag(1, :), ...
											handles.S.Fs, ...
											handles.CorrFrange, ...
											'Method', method, ...
											'Normalize', handles.NormalizeValue, ...
											'Lowcut', lowcut, ...
											'Level', handles.TargetSPL);
	end
	
	if strcmpi(handles.SignalMode, 'SYNTH')
		% apply ramp to stimuli
		handles.raw = sin2array(handles.raw, synth.ramp, handles.S.Fs);
		handles.adj = sin2array(handles.adj, synth.ramp, handles.S.Fs);
	end
	guidata(hObject, handles);
		
	% take fft of adj data
	[handles.fadj, handles.magadj, handles.phiadj] = ...
									daqdbfullfft(	handles.adj(bin(1):bin(2)), ...
														handles.S.Fs, ...
														length(handles.adj(bin(1):bin(2))) );
	guidata(hObject, handles);
	% update front panel plots
	updatePlots(hObject, handles);
	% save variables in workspace
	guidata(hObject, handles);
	disp('....updating complete');
%------------------------------------------------------------------------------

%------------------------------------------------------------------------------
function PlayRawSignalCtrl_Callback(hObject, eventdata, handles)
%------------------------------------------------------------------------------
	PlaySignal(hObject, handles);
%------------------------------------------------------------------------------

%------------------------------------------------------------------------------
function PlayAdjSignalCtrl_Callback(hObject, eventdata, handles)
%------------------------------------------------------------------------------
	PlaySignal(hObject, handles);
%------------------------------------------------------------------------------

%------------------------------------------------------------------------------
function LoadCalCtrl_Callback(hObject, eventdata, handles)
	% use the menu item callback
	LoadCalMenuItem_Callback(hObject, eventdata, handles)
%------------------------------------------------------------------------------

%------------------------------------------------------------------------------
function SaveSoundCtrl_Callback(hObject, eventdata, handles)
	% use the menu item callback
	SaveAdjSignalMenuItem_Callback(hObject, eventdata, handles);
%------------------------------------------------------------------------------

%------------------------------------------------------------------------------
function LoadWavCtrl_Callback(hObject, eventdata, handles)
	WavFilenameCtrl_Callback(hObject, eventdata, handles)
%------------------------------------------------------------------------------

%******************************************************************************
%******************************************************************************
%******************************************************************************

%******************************************************************************
%******************************************************************************
%******************************************************************************
% CALIBRATION DATA CONTROLS
%******************************************************************************
%******************************************************************************
%******************************************************************************
%------------------------------------------------------------------------------
function SmoothCalCtrl_Callback(hObject, eventdata, handles)
	% check value of handles.SmoothCalCtrl
	if read_ui_val(handles.SmoothCalCtrl)
		% if 1 (selected), enable the method text and control...
		enable_ui(handles.CalSmoothMethodText);
		enable_ui(handles.CalSmoothMethodCtrl);
		% ... and call the smoothmethod control callback
		CalSmoothMethodCtrl_Callback(hObject, eventdata, handles);
	else
		% if 0 (unselected), disable appropriate controls...
		disable_ui(handles.CalSmoothMethodText);
		disable_ui(handles.CalSmoothMethodCtrl);
		disable_ui(handles.SmoothVal1Text);
		disable_ui(handles.SmoothVal1Ctrl);
		% and plot original cal mag data
		plot(handles.CalibrationAxes, ...
				0.001*handles.cal.freq, ...
				handles.cal.mag(1, :), 'b.-');
		ylim( [	(0.9 * min(handles.cal.mag(1, :))) ...
					(1.1 * max(handles.cal.mag(1, :)))		]);
		grid on
	end
	guidata(hObject, handles);
%------------------------------------------------------------------------------

%------------------------------------------------------------------------------
function CalSmoothMethodCtrl_Callback(hObject, eventdata, handles)
%------------------------------------------------------------------------------
	% get value of smooth method
	smoothmethod = read_ui_val(handles.CalSmoothMethodCtrl);
	% based on smoothmethod, adjust UI controls/display
	switch(smoothmethod)
		case 1
			% moving-window average
			enable_ui(handles.SmoothVal1Text);
			update_ui_str(handles.SmoothVal1Text, 'Window Size');
			enable_ui(handles.SmoothVal1Ctrl);
			disable_ui(handles.SmoothVal2Text);
			disable_ui(handles.SmoothVal2Ctrl);
		case 2
			% Savitzky-Golay filter
			enable_ui(handles.SmoothVal1Text);
			update_ui_str(handles.SmoothVal1Text, 'Order');
			enable_ui(handles.SmoothVal1Ctrl);
			enable_ui(handles.SmoothVal2Text);
			update_ui_str(handles.SmoothVal2Text, 'Frame Size');
			enable_ui(handles.SmoothVal2Ctrl);
	end
	% smooth the calibration mag data, return in mag_smooth
	mag_smooth = SmoothCalibrationData(hObject, eventdata, handles);
	% plot the smoothed data
	plot(handles.CalibrationAxes, ...
				0.001*handles.cal.freq, ...
				mag_smooth(1, :), 'b.-');
	ylim([(0.9 * min(mag_smooth(1, :))) (1.1 * max(mag_smooth(1, :)))]);
	grid on
	% store smoothed data
	handles.mag_smooth = mag_smooth;
	guidata(hObject, handles);
%------------------------------------------------------------------------------

%------------------------------------------------------------------------------
function SmoothVal1Ctrl_Callback(hObject, eventdata, handles)
	tmp = read_ui_str(handles.SmoothVal1Ctrl, 'n');
	if tmp > 0
		handles.SmoothVal1 = tmp;
	else
		update_ui_str(handles.SmoothVal1Ctrl, handles.SmoothVal1);
		errordlg('Value must be greater than 0', 'Flat Wav Error');
	end
	CalSmoothMethodCtrl_Callback(hObject, eventdata, handles);
	guidata(hObject, handles)
%------------------------------------------------------------------------------

%------------------------------------------------------------------------------
function SmoothVal2Ctrl_Callback(hObject, eventdata, handles)
	tmp = read_ui_str(handles.SmoothVal2Ctrl, 'n');
	if tmp > 0
		if even(tmp)
			update_ui_str(handles.SmoothVal2Ctrl, handles.SmoothVal2);
			errordlg('Value must be odd', 'Flat Wav Error');
		else
			handles.SmoothVal2 = tmp;
		end
	else
		update_ui_str(handles.SmoothVal2Ctrl, handles.SmoothVal2);
		errordlg('Value must be greater than 0', 'Flat Wav Error');
	end
	CalSmoothMethodCtrl_Callback(hObject, eventdata, handles);
	guidata(hObject, handles)
%------------------------------------------------------------------------------

%------------------------------------------------------------------------------
function smoothed = SmoothCalibrationData(hObject, eventdata, handles)
	smoothmethod = read_ui_val(handles.CalSmoothMethodCtrl);
	switch(smoothmethod)
		case 1
			% moving window average
			[nrows, ncols] = size(handles.cal.mag);
			smoothed = zeros(nrows, ncols);
			% smooth each row of mags using moving_average() function
			% (internal to FlatWav)
			for n = 1:nrows
				smoothed(n, :) = moving_average(	handles.cal.mag(n, :), ...
															handles.SmoothVal1);
			end
		case 2
			% savitzky-golay filter
			[nrows, ncols] = size(handles.cal.mag);
			smoothed = zeros(nrows, ncols);
			for n = 1:nrows
				smoothed(n, :) = sgolayfilt(	handles.cal.mag(n, :), ...
														handles.SmoothVal1, ...
														handles.SmoothVal2	);
			end			
		
		otherwise
			% undefined method... should never get here if GUI is properly
			% functioning!!!
			smoothed = handles.cal.mag;
	end
%------------------------------------------------------------------------------

%------------------------------------------------------------------------------
function y = moving_average(x, w)
	%-----------------------------------------------------------------
	% computes sliding moving-average of vector x with window size w
	%-----------------------------------------------------------------
	k = ones(1, w) ./ w;
	% don't want to use normal 0 padding for x vector.  instead, add copies
	% of 1st and last values to x vector
	% first, force x to row vector
	x = x(:)';
	% store original length of x
	xlen = length(x);
	% then, pad vector with w copies of first and last x values
	x = [ x(1)*ones(1, w) x x(end)*ones(1, w)];
	% convolve
	y = conv(x, k, 'same');
	% truncate to proper length
	y = y(w + (1:xlen));
%------------------------------------------------------------------------------

%******************************************************************************
%******************************************************************************
%******************************************************************************


%******************************************************************************
%******************************************************************************
%******************************************************************************
% SELECT SIGNAL CONTROLS
%******************************************************************************
%******************************************************************************
%******************************************************************************

%------------------------------------------------------------------------------
function WavSignalButton_Callback(hObject, eventdata, handles)
	%--------------------------------------------------
	% user selected the Wav signal type
	%--------------------------------------------------
	% set signal mode to WAV
	handles.SignalMode = 'WAV';
	% make sure synth signal button is deselected
	update_ui_val(handles.SynthSignalButton, 0);
	% check if wavdata has been loaded/initialized
	if ~isempty(handles.wavdata.raw)
		% if so, copy raw signal and Fs to main handles 
		handles.raw = handles.wavdata.raw;
		handles.S.Fs = handles.wavdata.Fs;
	else
		% otherwise set handles.raw to empty
		handles.raw = [];
		handles.S.Fs = read_ui_str(handles.FsCtrl, 'n');
	end
	update_ui_val(handles.WavSignalButton, 1);
	guidata(hObject, handles);
	UpdateInputFilter(hObject, eventdata, handles);
%------------------------------------------------------------------------------

%------------------------------------------------------------------------------
function SynthSignalButton_Callback(hObject, eventdata, handles)
	%--------------------------------------------------
	% user selected the SYNTH signal type
	%--------------------------------------------------
	handles.SignalMode = 'SYNTH';
	update_ui_val(handles.WavSignalButton, 0);
	update_ui_val(handles.SynthSignalButton, 1);
	guidata(hObject, handles);
	updateGuiFromSynth(hObject, handles);
	UpdateInputFilter(hObject, eventdata, handles);
%------------------------------------------------------------------------------
%******************************************************************************
%******************************************************************************
%******************************************************************************


%******************************************************************************
%******************************************************************************
%******************************************************************************
% SYNTHESIZED SIGNAL CONTROLS
%******************************************************************************
%******************************************************************************
%******************************************************************************

%------------------------------------------------------------------------------
function SynthTypeCtrl_Callback(hObject, eventdata, handles)
	newVal = read_ui_val(handles.SynthTypeCtrl);
	% store old values
	switch handles.SynthIndex
		case 1
			handles.Tone = handles.synth;
		case 2
			handles.Noise = handles.synth;
		case 3
			handles.Sweep = handles.synth;
	end
	% update new ones
	switch newVal
		case 1
			handles.synth = handles.Tone;
		case 2
			handles.synth = handles.Noise;
		case 3
			handles.synth = handles.Sweep;
	end
	handles.SynthIndex = newVal;
	handles.SynthType = handles.S.Types{handles.SynthIndex};
	% store in handles
	guidata(hObject, handles);
	% update gui from synth settings
	updateGuiFromSynth(hObject, handles);
	UpdateInputFilter(hObject, eventdata, handles);
%------------------------------------------------------------------------------

%------------------------------------------------------------------------------
function FsCtrl_Callback(hObject, eventdata, handles)
%------------------------------------------------------------------------------
	Fs = read_ui_str(handles.FsCtrl, 'n');
	if between(Fs, 1, 1e6)
		handles.S.Fs = Fs;
	else
		update_ui_str(handles.FsCtrl, handles.S.Fs);
		disp('Bad sample rate')
	end
	guidata(hObject, handles);
	UpdateInputFilter(hObject, eventdata, handles);
%------------------------------------------------------------------------------

%------------------------------------------------------------------------------
function SynthCtrl_Callback(hObject, eventdata, handles)
%------------------------------------------------------------------------------
	% get the tag of the selected object
	tag = get(hObject, 'Tag');
	tagnum = find(strcmpi(tag, handles.S.CtrlTags));
	param = handles.S.Param{handles.SynthIndex}{tagnum};
	handles.synth.(param) = read_ui_str(hObject, 'n');
	guidata(hObject, handles);
	UpdateInputFilter(hObject, eventdata, handles);
%------------------------------------------------------------------------------
%******************************************************************************
%******************************************************************************
%******************************************************************************


%******************************************************************************
%******************************************************************************
%******************************************************************************
%% WAV FILE CONTROL
%******************************************************************************
%******************************************************************************
%******************************************************************************
%------------------------------------------------------------------------------
function WavFilenameCtrl_Callback(hObject, eventdata, handles)
%------------------------------------------------------------------------------
	[tmppath, tmpfile] = fileparts(handles.wavdata.datafile);
	[wavfile, wavpath] = uigetfile( '*.wav', ...
												'Load wav file...', ...
												tmppath);
	clear tmppath tmpfile
	if wavfile ~= 0
		wavdata.datafile = fullfile(wavpath, wavfile);
		[wavdata.raw, wavdata.Fs, wavdata.nbits, wavdata.opts] = wavread(wavdata.datafile);
		handles.wavdata = wavdata;
		update_ui_str(handles.WavFilenameCtrl, handles.wavdata.datafile);
		ostr = sprintf('Fs: %.2f\nnbits: %d\n', wavdata.Fs, wavdata.nbits);
		update_ui_str(handles.WaveInfoCtrl, ostr);
		handles.S.Fs = handles.wavdata.Fs;
		update_ui_str(handles.FsCtrl, handles.S.Fs);
		clear wavdata;
		% update input filter based on sample rate
		UpdateInputFilter(hObject, eventdata, handles);
		% update Analysis window based on duration of sound
		handles.Awindow(1) = 0;
		handles.Awindow(2) = floor(	bin2ms(length(handles.wavdata.raw), ...
												handles.S.Fs));
		update_ui_str(handles.AnalysisStartCtrl, handles.Awindow(1));
		update_ui_str(handles.AnalysisEndCtrl, handles.Awindow(2));
		guidata(hObject, handles);		
	else
		handles.wavdata = struct(	'datafile', '', 'raw', [], 'Fs', [], ...
											'nbits', [], 'opts', []);
		update_ui_str(handles.WavFilenameCtrl, '');
		update_ui_str(handles.WaveInfoCtrl, 'no wav loaded');
	end
	guidata(hObject, handles);
%-------------------------------------------------------------------------
%******************************************************************************
%******************************************************************************
%******************************************************************************


%******************************************************************************
%******************************************************************************
%******************************************************************************
% SIGNAL COMPENSATION CALLBACKS
%******************************************************************************
%******************************************************************************
%******************************************************************************

%------------------------------------------------------------------------------
% --- Executes on selection change in CompMethodCtrl.
function CompMethodCtrl_Callback(hObject, eventdata, handles)
	handles.CompMethod = read_ui_val(handles.CompMethodCtrl);
	guidata(hObject, handles);
%------------------------------------------------------------------------------

%------------------------------------------------------------------------------
function TargetSPLCtrl_Callback(hObject, eventdata, handles)
%------------------------------------------------------------------------------
	newVal = read_ui_str(handles.TargetSPLCtrl, 'n');
	% should perform some checks, forget for now
	if ~between(newVal, 0, 140)
		warndlg('Target SPL must be [0 140]');
		update_ui_str(handles.TargetSPL);
	else
		handles.TargetSPL = newVal;
	end
	guidata(hObject, handles);
%------------------------------------------------------------------------------

%------------------------------------------------------------------------------
function NormalizeCtrl_Callback(hObject, eventdata, handles)
%------------------------------------------------------------------------------
	newVal = read_ui_val(handles.NormalizeCtrl);
	if newVal
		handles.Normalize = 'on';
		enable_ui(handles.NormalizePeakCtrl);
		enable_ui(handles.NormalizePeakText);
	else
		handles.Normalize = 'off';
		disable_ui(handles.NormalizePeakCtrl);
		disable_ui(handles.NormalizePeakText);
	end
	guidata(hObject, handles);
%------------------------------------------------------------------------------

%------------------------------------------------------------------------------
function NormalizePeakCtrl_Callback(hObject, eventdata, handles)
%------------------------------------------------------------------------------
	newVal = read_ui_str(handles.NormalizePeakCtrl, 'n');
	if between(newVal, 0, 10)
		handles.NormalizeValue = newVal;
	else
		warndlg('Normalize value must be [0 10]');
		update_ui_str(handles.NormalizeValue);
	end
	guidata(hObject, handles);
%------------------------------------------------------------------------------

%------------------------------------------------------------------------------
function LowCutCtrl_Callback(hObject, eventdata, handles)
	newVal = read_ui_val(hObject);
	if newVal
		enable_ui(handles.LowCutFreqText);
		enable_ui(handles.LowCutFreqCtrl);
		handles.LowCut = 'on';
		handles.LowCutFreq = read_ui_str(handles.LowCutFreqCtrl, 'n');
	else
		disable_ui(handles.LowCutFreqText);
		disable_ui(handles.LowCutFreqCtrl);
		handles.LowCut = 'off';
	end
	guidata(hObject, handles);
%------------------------------------------------------------------------------

%------------------------------------------------------------------------------
function LowCutFreqCtrl_Callback(hObject, eventdata, handles)
	handles.LowCutFreq = read_ui_str(hObject, 'n');
	guidata(hObject, handles);
%------------------------------------------------------------------------------

%------------------------------------------------------------------------------
function CorrFminCtrl_Callback(hObject, eventdata, handles)
	newVal = read_ui_str(hObject, 'n');
	% NEED CHECKS!
	handles.CorrFrange(1) = newVal;
	guidata(hObject, handles);
%------------------------------------------------------------------------------

%------------------------------------------------------------------------------
function CorrFmaxCtrl_Callback(hObject, eventdata, handles)
	newVal = read_ui_str(hObject, 'n');
	% NEED CHECKS!
	handles.CorrFrange(2) = newVal;
	guidata(hObject, handles);
%------------------------------------------------------------------------------
%******************************************************************************
%******************************************************************************
%******************************************************************************

%******************************************************************************
%******************************************************************************
%******************************************************************************
% ANALYSIS WINDOW CONTROLS
%******************************************************************************
%******************************************************************************
%******************************************************************************

%------------------------------------------------------------------------------
function AnalysisStartCtrl_Callback(hObject, eventdata, handles)
	newVal = read_ui_str(hObject, 'n');
	% NEED CHECKS!
	handles.Awindow(1) = newVal;
	guidata(hObject, handles);
%------------------------------------------------------------------------------

%------------------------------------------------------------------------------
function AnalysisEndCtrl_Callback(hObject, eventdata, handles)
	newVal = read_ui_str(hObject, 'n');
	% NEED CHECKS!
	handles.Awindow(2) = newVal;
	guidata(hObject, handles);
%------------------------------------------------------------------------------

%******************************************************************************
%******************************************************************************
%******************************************************************************


%******************************************************************************
%******************************************************************************
%******************************************************************************
% OUTPUT DEVICE CONTROLS
%******************************************************************************
%******************************************************************************
%******************************************************************************

%------------------------------------------------------------------------------
% --- Executes on button press in SoundCardButton.
%------------------------------------------------------------------------------
function SoundCardButton_Callback(hObject, eventdata, handles)
	%--------------------------------------------------
	% user selected the SoundCard output
	%--------------------------------------------------
	% set signal output to WINSOUND
	handles.OutputDevice = 'WINSOUND';
	% make sure synth signal button is deselected
	update_ui_val(handles.NIDAQButton, 0);
	update_ui_val(handles.SoundCardButton, 1);
	update_ui_str(handles.RawdBText, sprintf('Raw dB SPL: ---'));
	update_ui_str(handles.AdjdBText, sprintf('Adj dB SPL: ---'));
	hide_uictrl(handles.RawdBText);
	hide_uictrl(handles.AdjdBText);
%	handles.S.Fs = read_ui_str(handles.FsCtrl, 'n');
	guidata(hObject, handles);
	UpdateInputFilter(hObject, eventdata, handles);
%------------------------------------------------------------------------------

%------------------------------------------------------------------------------
% --- Executes on button press in NIDAQButton.
%------------------------------------------------------------------------------
function NIDAQButton_Callback(hObject, eventdata, handles)
	%--------------------------------------------------
	% user selected the NIDAQ output
	%--------------------------------------------------
	% set signal output to NIDAQ
	handles.OutputDevice = 'NIDAQ';
	% make sure synth signal button is deselected
	update_ui_val(handles.NIDAQButton, 1);
	update_ui_val(handles.SoundCardButton, 0);
	show_uictrl(handles.RawdBText);
	show_uictrl(handles.AdjdBText);
	guidata(hObject, handles);
	UpdateInputFilter(hObject, eventdata, handles);
%------------------------------------------------------------------------------
%******************************************************************************
%******************************************************************************
%******************************************************************************


%******************************************************************************
%******************************************************************************
%******************************************************************************
%% Internal functions
%******************************************************************************
%******************************************************************************
%******************************************************************************

%------------------------------------------------------------------------------
function updateGuiFromSynth(hObject, handles)
	% get the index for the synthesized signal (ton = 1, noise = 2, sweep = 3) 
	sindx = handles.SynthIndex;
	% # of parameters for this signal type (for cycling through the S struct
	% values to be updated)
	Nsynthparam = handles.S.Nparam(sindx);
	% loop through parameters for this signal
	for n = 1:Nsynthparam
		update_ui_str(handles.(handles.S.TextTags{n}), [handles.S.Text{sindx}{n} ':']);
		update_ui_str(handles.(handles.S.CtrlTags{n}), handles.synth.(handles.S.Param{sindx}{n}));
		show_uictrl(handles.(handles.S.TextTags{n}));
		show_uictrl(handles.(handles.S.CtrlTags{n}));
	end
	% hide unused controls
	if handles.S.MaxNParam > Nsynthparam
		for n = (Nsynthparam+1):handles.S.MaxNParam
			hide_uictrl(handles.(handles.S.CtrlTags{n}));
			hide_uictrl(handles.(handles.S.TextTags{n}));
		end
	end
	% update ui elements
	update_ui_str(handles.FsCtrl, handles.S.Fs);
	update_ui_val(handles.SynthTypeCtrl, handles.SynthIndex);
	guidata(hObject, handles);
%------------------------------------------------------------------------------

%------------------------------------------------------------------------------
function updateSynthFromGui(hObject, handles)
	sindx = read_ui_val(handles.SynthTypeCtrl);
	handles.SynthIndex = sindx;
	handles.synth.type = handles.S.Types{sindx};
	for n = 1:handles.S.Nparam(sindx);
		handles.synth.(handles.S.Param{sindx}{n}) = ...
									read_ui_str(handles.(handles.S.CtrlTags{n}), 'n');
	end
	guidata(hObject, handles);
%------------------------------------------------------------------------------

%------------------------------------------------------------------------------
function updatePlots(hObject, handles)
	% plotting limits
	% limits for Mag and Phase plots
	dbmax = 0;
	dbmin = min([min(handles.magraw) min(handles.magadj)]);
	dblim = [dbmin dbmax];
	freqlim = 0.001*[0 handles.S.Fs/2];
	
	% update raw plots
	axes(handles.RawSignalAxes)
	tvec = 1000 * (0:(length(handles.raw)-1)) ./ handles.S.Fs;
	plot(tvec, handles.raw)
	title('Signal (V)')
	ylabel('Raw', 'Color', 'b')
	set(handles.RawSignalAxes, 'XTickLabel', []);
	xlim([min(tvec) max(tvec)])
	% get ticks
	time_ticks = get(handles.RawSignalAxes, 'XTick');
	
	axes(handles.RawMagAxes)
	plot(0.001*handles.fraw, handles.magraw);
	title('Magnitude (dB)')
	ylim(dblim);
	xlim(freqlim);
	set(handles.RawMagAxes, 'XTickLabel', []);
	
	axes(handles.RawPhaseAxes)
	plot(0.001*handles.fraw, unwrap(handles.phiraw));
	title('Phase (rad)')
	xlim(freqlim);
	set(handles.RawPhaseAxes, 'XTickLabel', []);
	
	axes(handles.RawSpectrumAxes)
	[S, F, T, P] = spectrogram(	handles.raw, ...
											handles.SpectrumWindow, ...
											[], ...
											handles.SpectrumWindow, ...
											handles.S.Fs	);
	save p.mat S F T P -MAT
	P = 20*log10(P);
	P(P == -Inf) = min(min(P(P ~= -Inf)));	
	surf(1000*T, 0.001*F, P, 'edgecolor', 'none');
	xlim([min(tvec) max(tvec)])
	ylim(freqlim);
	set(handles.RawSpectrumAxes, 'XTick', time_ticks)
	view(0, 90);
	title('Time vs. Freq (kHz) vs. dB')
	set(handles.RawSpectrumAxes, 'XTickLabel', []);
	colormap(handles.RawSpectrumAxes, handles.ColorMap)
% 	caxis([0.5*min(min(P)) max(max(P))])
	guidata(hObject, handles)
	
	% Update adj plots
	axes(handles.AdjSignalAxes)
	tvec = 1000 * (0:(length(handles.adj)-1)) ./ handles.S.Fs;
	plot(tvec, handles.adj, 'r')
	xlim([min(tvec) max(tvec)])
	ylabel('Adj', 'Color', 'r')
	xlabel('time (ms)')
	
	axes(handles.AdjMagAxes)
	plot(0.001*handles.fadj, handles.magadj, 'r');
	ylim(dblim);
	xlim(freqlim);
	xlabel('freq (kHz)');
	
	axes(handles.AdjPhaseAxes)
	plot(0.001*handles.fadj, unwrap(handles.phiadj), 'r');
	xlim(freqlim);
	xlabel('freq (kHz)');

	axes(handles.AdjSpectrumAxes)
% 	[S, F, T, P] = spectrogram(	handles.adj, ...
% 											handles.SpectrumWindow, ...
% 											floor(0.95*handles.SpectrumWindow), ...
% 											512, ...
% 											handles.S.Fs	);
	[S, F, T, P] = spectrogram(	handles.adj, ...
											handles.SpectrumWindow, ...
											[], ...
											handles.SpectrumWindow, ...
											handles.S.Fs	);
	P = 20*log10(P);
	P(P == -Inf) = min(min(P(P ~= -Inf)));	
	surf(1000*T, 0.001*F, P, 'edgecolor', 'none');
	xlim([min(tvec) max(tvec)])
	ylim(freqlim);
	set(handles.AdjSpectrumAxes, 'XTick', time_ticks)	
	view(0, 90);
	xlabel('Time (ms)')
	colormap(handles.AdjSpectrumAxes, handles.ColorMap)
% 	caxis([min(min(P)) max(max(P))])

	guidata(hObject, handles);
%------------------------------------------------------------------------------

%------------------------------------------------------------------------------
function [atten_val] = figure_atten(spl_val, rms_val, caldata)
%------------------------------------------------------------------------------
	[n, m] = size(caldata.mag);

	atten_val(1) = caldata.mindbspl(1) + db(rms_val(1)) - spl_val(1);
	if n == 2
		atten_val(2) = caldata.mindbspl(2) + db(rms_val(2)) - spl_val(2);
	end
%------------------------------------------------------------------------------


%******************************************************************************
%******************************************************************************
%******************************************************************************


%******************************************************************************
%******************************************************************************
%******************************************************************************
% MENU Callbacks
%******************************************************************************
%******************************************************************************
%******************************************************************************

%------------------------------------------------------------------------------
%------------------------------------------------------------------------------
% FlatWav Menu
%------------------------------------------------------------------------------
%------------------------------------------------------------------------------
%-------------------------------------------------------------------------
function SaveFigureMenuItem_Callback(hObject, eventdata, handles)
	[figfile, figpath] = uiputfile('*.fig','Save plot and figure in .fig file...');
	if figfile ~=0
		figfile = fullfile(figpath, figfile);
		saveas(handles.axes1, figfile, 'fig');
	end
%-------------------------------------------------------------------------

%-------------------------------------------------------------------------
function IndividualPlotMenuItem_Callback(hObject, eventdata, handles)
	% create new figure
	figure	
	% copy 
	a = axes;
	ax2ax(handles.axes1, a)
	plotStrings = read_ui_str(handles.PlotMenu);
	plotVal = read_ui_val(handles.PlotMenu);
	
	if isfield(handles.caldata, 'settings')
		if isfield(handles.caldata.settings, 'calfile')
			[pname, fname, tmp] = fileparts(handles.caldata.settings.calfile);
		else
			fname = {};
		end
	else
		fname = {};
	end	

	if ~isempty(fname)
		tstr = {plotStrings{plotVal}, fname};
	else
		tstr = plotStrings{plotVal};
	end
	
	title(tstr, 'Interpreter', 'none')
%-------------------------------------------------------------------------

%-------------------------------------------------------------------------
function PrintMenuItem_Callback(hObject, eventdata, handles)
%-------------------------------------------------------------------------
	printdlg(handles.figure1)
%-------------------------------------------------------------------------

%-------------------------------------------------------------------------
function QuitMenuItem_Callback(hObject, eventdata, handles)
%-------------------------------------------------------------------------
	selection = questdlg(['Close ' get(handles.figure1,'Name') '?'],...
								['Close ' get(handles.figure1,'Name') '...'],...
								'Yes','No','Yes');
	if strcmp(selection,'No')
		 return;
	end
	delete(handles.figure1);
%------------------------------------------------------------------------------

%------------------------------------------------------------------------------
%------------------------------------------------------------------------------
% Cal Menu
%------------------------------------------------------------------------------
%------------------------------------------------------------------------------

%-------------------------------------------------------------------------
function LoadCalMenuItem_Callback(hObject, eventdata, handles)
%-------------------------------------------------------------------------
	[calfile, calpath] = uigetfile( {'*.cal'; '*_cal.mat'}, ...
												'Load headphone calibration data from file...');
	if calfile ~=0
		datafile = fullfile(calpath, calfile);	
		handles.cal = load_headphone_cal(datafile);
		plot(	handles.CalibrationAxes, ...
				0.001*handles.cal.freq, ...
				handles.cal.mag(1, :), '.-');
		ylim([0.9*min(handles.cal.mag(1, :)) 1.1*max(handles.cal.mag(1, :))]);
		grid on
	end
	guidata(hObject, handles);
%-------------------------------------------------------------------------


%-------------------------------------------------------------------------
function FlatCalMenuItem_Callback(hObject, eventdata, handles)
%--------------------------------------------------
% fake cal data
%--------------------------------------------------
	handles.cal = fake_caldata('freqs', (1:10:(handles.S.Fs / 2)));
	handles.cal.mag = 90 * handles.cal.mag;
	guidata(hObject, handles);
	plot(handles.CalibrationAxes, 0.001*handles.cal.freq, handles.cal.mag(1, :), '.-');
	ylim([0 100]);
%-------------------------------------------------------------------------

%------------------------------------------------------------------------------
%------------------------------------------------------------------------------
% Wav Menu
%------------------------------------------------------------------------------
%------------------------------------------------------------------------------

%-------------------------------------------------------------------------
function LoadWavMenuItem_Callback(hObject, eventdata, handles)
%-------------------------------------------------------------------------
	WavFilenameCtrl_Callback(hObject, eventdata, handles);
%-------------------------------------------------------------------------	

%------------------------------------------------------------------------------
%------------------------------------------------------------------------------
% Signal Menu
%------------------------------------------------------------------------------
%------------------------------------------------------------------------------

%-------------------------------------------------------------------------
function SaveAdjSignalMenuItem_Callback(hObject, eventdata, handles)
%-------------------------------------------------------------------------
	if isempty(handles.adj)
		warndlg(sprintf('%s: adj vector is empty!  aborting save', mfilename));
		return
	end
	
	[tmppath, tmpname, tmpext] = fileparts(handles.wavdata.datafile);
	tmpfile = fullfile(tmppath, [tmpname '_adj' tmpext]);
	
	[adjfile, adjpath] = uiputfile(	'*.wav', ...
												'Save adj signal to wav file...', tmpfile);
	if adjfile ~=0
		datafile = fullfile(adjpath, adjfile);
		peakval = max(handles.adj);
		if peakval >= 1
			fprintf('!!!!!!!\nPoints in adj are >= 1\nFile will be normalized\n');
			wavwrite(0.9*normalize(handles.adj), handles.S.Fs, datafile);
% 			peakfile = [datafile(1:(end-4)) '_PeakVal.txt'];
% 			save(peakfile, peakval, '-ascii');
		else
			wavwrite(handles.adj, handles.S.Fs, datafile);
		end
	end
	guidata(hObject, handles);
%-------------------------------------------------------------------------

%-------------------------------------------------------------------------
function SaveRawSignalMenuItem_Callback(hObject, eventdata, handles)
%-------------------------------------------------------------------------
	if isempty(handles.raw)
		warndlg(sprintf('%s: raw vector is empty!  aborting save', mfilename));
		return
	end
	
	[rawfile, rawpath] = uiputfile(	'*.wav', ...
												'Save raw signal to wav file...');
	if rawfile ~=0
		datafile = fullfile(rawpath, rawfile);
		peakval = max(handles.raw);
		if peakval >= 1
			fprintf('!!!!!!!\nPoints in raw are >= 1\nFile will be normalized\n');
			wavwrite(0.9*normalize(handles.raw), handles.S.Fs, datafile);
			peakfile = [datafile(1:(end-4)) '_PeakVal.txt'];
			save(peakfile, peakval, '-ascii');
		else
			wavwrite(handles.raw, handles.S.Fs, datafile);
		end
	end
%-------------------------------------------------------------------------

%-------------------------------------------------------------------------
function SaveAllSignalsMenuItem_Callback(hObject, eventdata, handles)
%-------------------------------------------------------------------------
	[matfile, matpath] = uiputfile(	'*.mat', ...
												'Save signals to mat file...');
	if matfile ~= 0
		raw = handles.raw;
		adj = handles.adj;
		S = handles.S;
		if strcmpi(handles.SignalMode, 'WAV')
			wavdata = handles.wavdata;
			save(fullfile(matpath, matfile), 'raw', 'adj', 'S', 'wavdata', '-MAT');
			clear raw adj S wavdata;
		else
			synth = handles.synth;
			save(fullfile(matpath, matfile), 'raw', 'adj', 'Fs', 'synth', '-MAT');
			clear raw adj S synth;
		end
	end
%-------------------------------------------------------------------------

%------------------------------------------------------------------------------
%------------------------------------------------------------------------------
% SETTINGS Menu
%------------------------------------------------------------------------------
%------------------------------------------------------------------------------

%-------------------------------------------------------------------------
%-------------------------------------------------------------------------
function MicSensitivityMenuItem_Callback(hObject, eventdata, handles)
	newVal = uiaskvalue(	'Value',				handles.MicSensitivity,			...
								'ValueText',		'Mic Sensitivity (V/Pa)',		...
								'QuestionText',	'Enter value from NEXXUS',		...
								'FigureName',		''	);
	handles.MicSensitivity = newVal;
	% pre-compute the V -> Pa conversion factor
	handles.VtoPa = (1/handles.MicGain) * (1/handles.MicSensitivity);
	guidata(hObject, handles);
%-------------------------------------------------------------------------
%-------------------------------------------------------------------------

%-------------------------------------------------------------------------
%-------------------------------------------------------------------------
function MicrophoneGainMenuItem_Callback(hObject, eventdata, handles)
	newVal = uiaskvalue(	'Value',				handles.MicGaindB,			...
								'ValueText',		'Mic Gain (dB)',		...
								'QuestionText',	'Enter mic gain (0 dB for NEXXUS)',...
								'FigureName',		''	);
	handles.MicGaindB = newVal;
	handles.MicGain = invdb(handles.MicGaindB);
	% pre-compute the V -> Pa conversion factor
	handles.VtoPa = (1/handles.MicGain) * (1/handles.MicSensitivity);
	guidata(hObject, handles);
%-------------------------------------------------------------------------
%-------------------------------------------------------------------------

%-------------------------------------------------------------------------
%-------------------------------------------------------------------------
function InputFilterMenuItem_Callback(hObject, eventdata, handles)
	% get new low pass cutoff frequency
	newVal = uiaskvalue(	...
					'Value',				handles.LPFc,			...
					'ValueText',		'Lowpass Filter Fc (Hz)',		...
					'QuestionText',	'Input LowPass Filter Cutoff Frequency', ...
					'FigureName',		''	);
	handles.LPFc = newVal;
	% get new highpass cutoff frequency
	newVal = uiaskvalue(	...
					'Value',				handles.HPFc,			...
					'ValueText',		'Highpass Filter Fc (Hz)',		...
					'QuestionText',	'Input HighPass Filter Cutoff Frequency', ...
					'FigureName',		''	);
	%{
	handles.HPFc = newVal;
	newVal = uiaskvalue(	'Value',				handles.FilterOrder,			...
								'ValueText',		'Filter Order (>0)',		...
								'QuestionText',	'Input Filter Order', ...
								'FigureName',		''	);
	handles.FilterOrder = newVal;
	%}

	%--------------------------------------------------------------
	% Define new bandpass filter coeffs for processing the data
	%--------------------------------------------------------------
	% Nyquist frequency
	fnyq = handles.S.Fs / 2;
	% passband definition
	fband = [handles.HPFc handles.LPFc] ./ fnyq;
	% filter coefficients using a butterworth bandpass filter
	[handles.fcoeffb, handles.fcoeffa] = butter( handles.FilterOrder, ...
																fband, 'bandpass');
	% store settings
	guidata(hObject, handles);
%-------------------------------------------------------------------------

%-------------------------------------------------------------------------
%-------------------------------------------------------------------------
function SpectrumWindowMenuItem_Callback(hObject, eventdata, handles)
	% get new spectrum window
	newVal = uiaskvalue(	...
					'Value',				handles.SpectrumWindow,			...
					'ValueText',		'Spectrum Window (pts)',		...
					'QuestionText',	'Input Window Size for Spectrum', ...
					'FigureName',		''	);
	if ~isnumeric(newVal)
		warndlg('Spectrum Window size must be a number!', 'FlatWav')
	elseif ~between(newVal, 2, 1e6)
		warndlg('Spectrum Window size must be between 2 and 1e6!', 'FlatWav')
	else
		handles.SpectrumWindow = newVal;
	end
	guidata(hObject, handles);
	updatePlots(hObject, handles);
%-------------------------------------------------------------------------

%******************************************************************************
%******************************************************************************
%******************************************************************************

%******************************************************************************
%******************************************************************************
%******************************************************************************
% Executes during object creation, after setting all properties.
%******************************************************************************
%******************************************************************************
%******************************************************************************
function CompMethodCtrl_CreateFcn(hObject, eventdata, handles)
	if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
		  set(hObject,'BackgroundColor','white');
	end
	set(hObject, 'String', {'plot(rand(5))', 'plot(sin(1:0.01:25))', 'bar(1:.5:10)', 'plot(membrane)', 'surf(peaks)'});
function WavFilenameCtrl_CreateFcn(hObject, eventdata, handles)
	if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
		 set(hObject,'BackgroundColor','white');
	end
function SynthTypeCtrl_CreateFcn(hObject, eventdata, handles)
	if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
		 set(hObject,'BackgroundColor','white');
	end
function FsCtrl_CreateFcn(hObject, eventdata, handles)
	if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
		 set(hObject,'BackgroundColor','white');
	end
function p1Ctrl_CreateFcn(hObject, eventdata, handles)
	if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
		 set(hObject,'BackgroundColor','white');
	end
function p2Ctrl_CreateFcn(hObject, eventdata, handles)
	if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
		 set(hObject,'BackgroundColor','white');
	end
function p3Ctrl_CreateFcn(hObject, eventdata, handles)
	if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
		 set(hObject,'BackgroundColor','white');
	end
function p4Ctrl_CreateFcn(hObject, eventdata, handles)
	if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
		 set(hObject,'BackgroundColor','white');
	end
function p5Ctrl_CreateFcn(hObject, eventdata, handles)
	if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
		 set(hObject,'BackgroundColor','white');
	end
function p6Ctrl_CreateFcn(hObject, eventdata, handles)
	if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
		 set(hObject,'BackgroundColor','white');
	end
function TargetSPLCtrl_CreateFcn(hObject, eventdata, handles)
	if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
		 set(hObject,'BackgroundColor','white');
	end
function LowCutFreqCtrl_CreateFcn(hObject, eventdata, handles)
	if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
		 set(hObject,'BackgroundColor','white');
	end
function SpectrumWindowCtrl_CreateFcn(hObject, eventdata, handles)
	if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
		 set(hObject,'BackgroundColor','white');
	end
function CorrFminCtrl_CreateFcn(hObject, eventdata, handles)
	if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
		 set(hObject,'BackgroundColor','white');
	end
function CorrFmaxCtrl_CreateFcn(hObject, eventdata, handles)
	if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
		 set(hObject,'BackgroundColor','white');
	end
function NormalizePeakCtrl_CreateFcn(hObject, eventdata, handles)
	if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
		 set(hObject,'BackgroundColor','white');
	end
function AnalysisStartCtrl_CreateFcn(hObject, eventdata, handles)
	if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
		 set(hObject,'BackgroundColor','white');
	end
function AnalysisEndCtrl_CreateFcn(hObject, eventdata, handles)
	if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
		 set(hObject,'BackgroundColor','white');
	end
function SmoothVal1Ctrl_CreateFcn(hObject, eventdata, handles)
	if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
		 set(hObject,'BackgroundColor','white');
	end

function CalSmoothMethodCtrl_CreateFcn(hObject, eventdata, handles)
	if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
		set(hObject,'BackgroundColor','white');
	end
function SmoothVal2Ctrl_CreateFcn(hObject, eventdata, handles)
	if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
	    set(hObject,'BackgroundColor','white');
	end
%******************************************************************************
%******************************************************************************
%******************************************************************************








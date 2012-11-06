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

% Last Modified by GUIDE v2.5 22-Oct-2012 18:37:14

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

%------------------------------------------------------------------------------
%------------------------------------------------------------------------------
%------------------------------------------------------------------------------
%% Essential Functions
%------------------------------------------------------------------------------
%------------------------------------------------------------------------------
%------------------------------------------------------------------------------

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
	% directory when using installed version:
	%	pdir = ['C:\TytoLogy\TytoLogySettings\' getenv('USERNAME')];
	% development tree
	pdir = ['C:\Users\sshanbhag\Code\Matlab\TytoLogy\TytoLogySettings\' getenv('USERNAME')];
	if isempty(which('ms2samples'))
		% could not find the RPload.m function (which is in TytoLogy
		% toolbox) which suggests that the paths are not set or are 
		% incorrect for this setup.  load the paths using the tytopaths program.
		%--------
		% First, store the current path
		cdir = pwd;
		% build the path to the user's TytoSettings directory and
		% change dirs to it.  Run the tytopaths script and then
		% return to the original ("current") directory
		disp([mfilename ': loading paths using ' pdir])
		cd(pdir);
		tytopaths
		cd(cdir);
	else
		disp([mfilename ': paths ok, launching programn'])
	end

	%--------------------------------------------------
	% SYNTH SETTINGS
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

	%--------------------------------------------------
	% update GUI and synth
	%--------------------------------------------------
	guidata(hObject, handles);
	updateGuiFromSynth(hObject, handles)
	guidata(hObject, handles);
	updateSynthFromGui(hObject, handles);
	guidata(hObject, handles);

	%--------------------------------------------------
	% COMPENSATION SETTINGS
	%--------------------------------------------------
	% reset string in CompMethodCtrl
	set(handles.CompMethodCtrl, 'string', 'none|atten|boost|compress');
	% set compensation method to 1 ('atten')
	handles.CompMethod = 1;
	update_ui_val(handles.CompMethodCtrl, handles.CompMethod);
	% set correction freq range (in Hz)
	handles.CorrFrange = [10 10000];
	% default normalize status
	handles.Normalize = 'on';
	handles.NormalizeValue = 1.0;
	update_ui_val(handles.NormalizeCtrl, handles.NormalizeValue);
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
	% spectrum settings
	%--------------------------------------------------
	handles.SpectrumWindow = 1024;
	update_ui_str(handles.SpectrumWindowCtrl, handles.SpectrumWindow);
	guidata(hObject, handles);

	%--------------------------------------------------
	% set initial state for sounds
	%--------------------------------------------------
	% wavdata struct holds information about wav file.
	% create blank wavdata struct, update gui
	handles.wavdata = struct(	'datafile', [], 'raw', [], 'fs', [], ...
										'nbits', [], 'opts', []);
	update_ui_str(handles.FilenameCtrl, '');
	update_ui_str(handles.WaveInfoCtrl, 'no wav loaded');
	% create empty raw and adj vectors
	handles.raw = [];
	handles.adj = [];
	guidata(hObject, handles);

	%--------------------------------------------------
	% fake cal data
	%--------------------------------------------------
	handles.cal = fake_caldata('freqs', [1:10:(handles.S.Fs / 2)]);
	handles.cal.mag = 90 * handles.cal.mag;
	guidata(hObject, handles);
	plot(handles.CalibrationAxes, 0.001*handles.cal.freq, handles.cal.mag(1, :), '.-');
	ylim([0 100]);

%------------------------------------------------------------------------------
%------------------------------------------------------------------------------

%------------------------------------------------------------------------------
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
%------------------------------------------------------------------------------

%------------------------------------------------------------------------------
%------------------------------------------------------------------------------
%------------------------------------------------------------------------------
%% Ctrl Callbacks
%------------------------------------------------------------------------------
%------------------------------------------------------------------------------
%------------------------------------------------------------------------------

%------------------------------------------------------------------------------
% --- Executes on selection change in CompMethodCtrl.
function CompMethodCtrl_Callback(hObject, eventdata, handles)
	handles.CompMethod = read_ui_val(handles.CompMethodCtrl);
	guidata(hObject, handles);
%------------------------------------------------------------------------------


%------------------------------------------------------------------------------
% --- Executes on button press in UpdateSignalCtrl.
%------------------------------------------------------------------------------
function UpdateSignalCtrl_Callback(hObject, eventdata, handles)
%------------------------------------------------------------------------------
% using settings from GUI, update the test signal, compensated signal, and
% plots of the signals
%------------------------------------------------------------------------------
	%--------------------------------------------------------------------
	% Signal Mode is either WAV (.wav file) or SYNTH (synthesized)
	%--------------------------------------------------------------------
	switch handles.SignalMode
		
		case 'WAV'
			% load wav file
			[wavdata.raw, wavdata.Fs, wavdata.nbits, wavdata.opts] = wavread(handles.wavdata.datafile);
			% apply ramp (short, just to ensure zeros at beginning and end of
			% stimulus)
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
					handles.raw = synmonosine(synth.dur, handles.S.Fs, synth.freq, synth.amp, 0);

				case 'noise'
					% create noise
					handles.raw = synmononoise_fft(synth.dur, handles.S.Fs, synth.fmin, synth.fmax, synth.amp, 0);
					% kludge to scale amplitude properly
					handles.raw = synth.amp * normalize(handles.raw);
					
				case 'sweep'
					% create FM sweep (via wrapper around chirp() matlab
					% function)
					handles.raw = synsweep(synth.dur, handles.S.Fs, synth.fmin, synth.fmax, synth.amp, 0);
			end
			% apply ramp to raw stimulus
			handles.raw = sin2array(handles.raw, synth.ramp, handles.S.Fs);
	end
	
	% take fft of raw data
	[handles.fraw, handles.magraw, handles.phiraw] = ...
											daqdbfullfft(handles.raw, handles.S.Fs, length(handles.raw));
	guidata(hObject, handles);

	%--------------------------------------------------------------------
	% apply compensation method
	% CompMethod is value of 1, 2, 3, 4 (value of CompMethodCtrl) which 
	% is a pull-down menu.  Map these to string values
	%--------------------------------------------------------------------
	switch handles.CompMethod
		case 1
			method = 'NONE'
		case 2
			method = 'ATTEN'
		case 3
			method = 'BOOST'
		case 4
			method = 'COMPRESS'
	end
	
	if strcmp(handles.LowCut, 'off')
		lowcut = 'off';
	else
		lowcut = handles.LowCutFreq;
	end
	
	if strcmpi(method, 'NONE')
		handles.adj = handles.raw;
		
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
	
	% take fft of adj data
	[handles.fadj, handles.magadj, handles.phiadj] = daqdbfullfft(handles.adj, handles.S.Fs, length(handles.adj));
	guidata(hObject, handles);

	updatePlots(hObject, handles);
	guidata(hObject, handles);
%------------------------------------------------------------------------------

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
%------------------------------------------------------------------------------

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

	guidata(hObject, handles);
	updateGuiFromSynth(hObject, handles);
		
%------------------------------------------------------------------------------

%------------------------------------------------------------------------------
function FsCtrl_Callback(hObject, eventdata, handles)
	Fs = read_ui_str(handles.FsCtrl, 'n');
	if between(Fs, 1, 1e6)
		handles.S.Fs = Fs;
	else
		update_ui_str(handles.FsCtrl, handles.S.Fs);
		disp('Bad sample rate')
	end
	guidata(hObject, handles);
%------------------------------------------------------------------------------

%------------------------------------------------------------------------------
function SynthCtrl_Callback(hObject, eventdata, handles)
	% get the tag of the selected object
	tag = get(hObject, 'Tag');
	tagnum = find(strcmpi(tag, handles.S.CtrlTags));
	param = handles.S.Param{handles.SynthIndex}{tagnum}
	handles.synth.(param) = read_ui_str(hObject, 'n');
	guidata(hObject, handles);
%------------------------------------------------------------------------------

%------------------------------------------------------------------------------
function FilenameCtrl_Callback(hObject, eventdata, handles)
	loadWavFile(hObject, eventdata, handles);
%-------------------------------------------------------------------------

%------------------------------------------------------------------------------

%------------------------------------------------------------------------------
function TargetSPLCtrl_Callback(hObject, eventdata, handles)
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
function PlaySignalCtrl_Callback(hObject, eventdata, handles)
	if (handles.S.Fs == 44100)  && ~isempty(handles.raw) && ~isempty(handles.adj)
		disable_ui(hObject);
		update_ui_str(hObject, 'RAW');
		pause(0.5);
		sound(handles.raw, handles.S.Fs);
		pause(2);
		update_ui_str(hObject, 'ADJ');
		pause(0.5);
		sound(handles.adj, handles.S.Fs);
		pause(1);
		enable_ui(hObject);
		update_ui_str(hObject, 'Play');
	end
%------------------------------------------------------------------------------

%------------------------------------------------------------------------------
function NormalizeCtrl_Callback(hObject, eventdata, handles)
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
function SaveSoundCtrl_Callback(hObject, eventdata, handles)
	% use the menu item callback
	SaveAdjSignalMenuItem_Callback(hObject, eventdata, handles);
%------------------------------------------------------------------------------

%------------------------------------------------------------------------------
function LoadCalCtrl_Callback(hObject, eventdata, handles)
	% use the menu item callback
	LoadCalMenuItem_Callback(hObject, eventdata, handles)
%------------------------------------------------------------------------------

%------------------------------------------------------------------------------
function SpectrumWindowCtrl_Callback(hObject, eventdata, handles)
	newVal = read_ui_str(handles.SpectrumWindowCtrl, 'n');
	if ~isnumeric(newVal)
		update_ui_str(handles.SpectrumWindowCtrl, handles.SpectrumWindow);
		warning('Spectrum Window size must be a number!')
	elseif ~between(newVal, 2, 1e6)
		update_ui_str(handles.SpectrumWindowCtrl, handles.SpectrumWindow);
		warning('Spectrum Window size must be between 2 and 1e6!')
	else
		handles.SpectrumWindow = newVal;
	end
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

%------------------------------------------------------------------------------
%------------------------------------------------------------------------------
%------------------------------------------------------------------------------
%% Internal functions
%------------------------------------------------------------------------------
%------------------------------------------------------------------------------
%------------------------------------------------------------------------------

%------------------------------------------------------------------------------
function updateGuiFromSynth(hObject, handles)
	sindx = handles.SynthIndex;
	Nsynthparam = handles.S.Nparam(sindx);
	
	for n = 1:Nsynthparam
		update_ui_str(handles.(handles.S.TextTags{n}), [handles.S.Text{sindx}{n} ':']);
		update_ui_str(handles.(handles.S.CtrlTags{n}), handles.synth.(handles.S.Param{sindx}{n}));
		show_uictrl(handles.(handles.S.TextTags{n}));
		show_uictrl(handles.(handles.S.CtrlTags{n}));
	end
	if handles.S.MaxNParam > Nsynthparam
		for n = (Nsynthparam+1):handles.S.MaxNParam
			hide_uictrl(handles.(handles.S.CtrlTags{n}));
			hide_uictrl(handles.(handles.S.TextTags{n}));
		end
	end
	update_ui_str(handles.FsCtrl, handles.S.Fs);
	update_ui_val(handles.SynthTypeCtrl, handles.SynthIndex);
	
%------------------------------------------------------------------------------

%------------------------------------------------------------------------------
function updateSynthFromGui(hObject, handles)
	sindx = read_ui_val(handles.SynthTypeCtrl);
	handles.SynthIndex = sindx;
	handles.synth.type = handles.S.Types{sindx};
	for n = 1:handles.S.Nparam(sindx);
		handles.synth.(handles.S.Param{sindx}{n}) = read_ui_str(handles.(handles.S.CtrlTags{n}), 'n');
	end
	guidata(hObject, handles);
%------------------------------------------------------------------------------

function updatePlots(hObject, handles)
	% plotting limits
	dblim = [-120 0];
	freqlim = 0.001*[0 handles.S.Fs/2];

	% update raw plots
	axes(handles.RawSignalAxes)
	tvec = 1000 * (0:(length(handles.raw)-1)) ./ handles.S.Fs;
	plot(tvec, handles.raw)
	title('Signal (V)')
	set(handles.RawSignalAxes, 'XTickLabel', []);
	
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
											floor(0.95*handles.SpectrumWindow), ...
											512, ...
											handles.S.Fs	);
	surf(1000*T, 0.001*F, 20*log10(P), 'edgecolor', 'none');
	ylim(freqlim);
	axis tight;
	view(0, 90);
	title('Time vs. Freq (kHz) vs. dB')
	set(handles.RawSpectrumAxes, 'XTickLabel', []);
	colormap('gray')
	
	% Update adj plots
	axes(handles.AdjSignalAxes)
	tvec = 1000 * (0:(length(handles.adj)-1)) ./ handles.S.Fs;
	plot(tvec, handles.adj)
	xlabel('time (ms)')
	
	axes(handles.AdjMagAxes)
	plot(0.001*handles.fadj, handles.magadj);
	ylim(dblim);
	xlim(freqlim);
	xlabel('freq (kHz)');
	
	axes(handles.AdjPhaseAxes)
	plot(0.001*handles.fadj, unwrap(handles.phiadj));
	xlim(freqlim);
	xlabel('freq (kHz)');

	axes(handles.AdjSpectrumAxes)
	[S, F, T, P] = spectrogram(	handles.adj, ...
											handles.SpectrumWindow, ...
											floor(0.95*handles.SpectrumWindow), ...
											512, ...
											handles.S.Fs	);
	surf(1000*T, 0.001*F, 20*log10(P), 'edgecolor', 'none');
	ylim(freqlim);
	axis tight;
	view(0, 90);
	xlabel('Time (ms)')
	colormap('gray')
	
	
	guidata(hObject, handles);
%------------------------------------------------------------------------------

%-------------------------------------------------------------------------
function loadWavFile(hObject, eventdata, handles)
	[wavfile, wavpath] = uigetfile( '*.wav', ...
												'Load wav file...');
	if wavfile ~= 0
		wavdata.datafile = fullfile(wavpath, wavfile);
		[wavdata.raw, wavdata.fs, wavdata.nbits, wavdata.opts] = wavread(wavdata.datafile);
		handles.wavdata = wavdata;
		update_ui_str(handles.FilenameCtrl, wavdata.datafile);
		update_ui_str(handles.WaveInfoCtrl, wavdata.datafile);
		clear wavdata;
		
	else
		handles.wavdata = struct(	'datafile', [], 'raw', [], 'fs', [], ...
											'nbits', [], 'opts', []);
		update_ui_str(handles.FilenameCtrl, '');
		update_ui_str(handles.WaveInfoCtrl, 'no wav loaded');
	end
	guidata(hObject, handles);
%-------------------------------------------------------------------------

%------------------------------------------------------------------------------
%------------------------------------------------------------------------------
%------------------------------------------------------------------------------
%% MENU Callbacks
%------------------------------------------------------------------------------
%------------------------------------------------------------------------------
%------------------------------------------------------------------------------

%------------------------------------------------------------------------------
function CloseMenuItem_Callback(hObject, eventdata, handles)
	selection = questdlg(['Close ' get(handles.figure1,'Name') '?'],...
								['Close ' get(handles.figure1,'Name') '...'],...
								'Yes','No','Yes');
	if strcmp(selection,'No')
		 return;
	end

	delete(handles.figure1)
%------------------------------------------------------------------------------

%-------------------------------------------------------------------------
function LoadCalMenuItem_Callback(hObject, eventdata, handles)
	[calfile, calpath] = uigetfile( {'*.cal'; '*_cal.mat'}, ...
												'Load headphone calibration data from file...');
	if calfile ~=0
		datafile = fullfile(calpath, calfile);	
		handles.cal = load_headphone_cal(datafile);
		plot(handles.CalibrationAxes, 0.001*handles.cal.freq, handles.cal.mag(1, :), '.-');
		ylim([min(handles.cal.mag(1, :)) max(handles.cal.mag(1, :))]);
		grid on
	end
	guidata(hObject, handles);
%-------------------------------------------------------------------------

%-------------------------------------------------------------------------
function PrintMenuItem_Callback(hObject, eventdata, handles)
	printdlg(handles.figure1)
%-------------------------------------------------------------------------

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
			[pname, fname, ~] = fileparts(handles.caldata.settings.calfile);
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
function LoadWavMenuItem_Callback(hObject, eventdata, handles)
	loadWavFile(hObject, eventdata, handles);

%-------------------------------------------------------------------------

%-------------------------------------------------------------------------
function SaveAdjSignalMenuItem_Callback(hObject, eventdata, handles)
	if isempty(handles.adj)
		warning('%s: adj vector is empty!  aborting save', mfilename);
		return
	end
	
	[adjfile, adjpath] = uiputfile(	'*.wav', ...
												'Save adj signal to wav file...');
	if adjfile ~=0
		datafile = fullfile(adjpath, adjfile);
		peakval = max(handles.adj);
		if peakval >= 1
			fprintf('!!!!!!!!!!!!!!!!\nPoints in adj are >= 1\nFile will be normalized\n');
			wavwrite(0.9*normalize(handles.adj), handles.S.Fs, datafile);
			peakfile = [datafile(1:(end-4)) '_PeakVal.txt'];
			save(peakfile, peakval, '-ascii');
		else
			wavwrite(handles.adj, handles.S.Fs, datafile);
		end
	end
%-------------------------------------------------------------------------

%-------------------------------------------------------------------------
function SaveRawSignalMenuItem_Callback(hObject, eventdata, handles)
	if isempty(handles.raw)
		warning('%s: raw vector is empty!  aborting save', mfilename);
		return
	end
	
	[rawfile, rawpath] = uiputfile(	'*.wav', ...
												'Save raw signal to wav file...');
	if rawfile ~=0
		datafile = fullfile(rawpath, rawfile);
		peakval = max(handles.raw);
		if peakval >= 1
			fprintf('!!!!!!!!!!!!!!!!\nPoints in raw are >= 1\nFile will be normalized\n');
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

%-------------------------------------------------------------------------
function FlatCalMenuItem_Callback(hObject, eventdata, handles)
%--------------------------------------------------
% fake cal data
%--------------------------------------------------
	handles.cal = fake_caldata('freqs', [1:10:(handles.S.Fs / 2)]);
	handles.cal.mag = 90 * handles.cal.mag;
	guidata(hObject, handles);
	plot(handles.CalibrationAxes, 0.001*handles.cal.freq, handles.cal.mag(1, :), '.-');
	ylim([0 100]);
%-------------------------------------------------------------------------

%------------------------------------------------------------------------------
%------------------------------------------------------------------------------
%------------------------------------------------------------------------------


%------------------------------------------------------------------------------
%------------------------------------------------------------------------------
%------------------------------------------------------------------------------
%% Executes during object creation, after setting all properties.
%------------------------------------------------------------------------------
%------------------------------------------------------------------------------
function CompMethodCtrl_CreateFcn(hObject, eventdata, handles)
	if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
		  set(hObject,'BackgroundColor','white');
	end
	set(hObject, 'String', {'plot(rand(5))', 'plot(sin(1:0.01:25))', 'bar(1:.5:10)', 'plot(membrane)', 'surf(peaks)'});
function FilenameCtrl_CreateFcn(hObject, eventdata, handles)
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
%------------------------------------------------------------------------------
%------------------------------------------------------------------------------
%------------------------------------------------------------------------------





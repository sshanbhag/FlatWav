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

% Last Modified by GUIDE v2.5 12-Oct-2012 12:33:56

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
% reset string
set(handles.CompMethodCtrl, 'string', 'none|atten|boost|compress');
% set compensation method to 1 ('atten')
handles.CompMethod = 1;
update_ui_val(handles.CompMethodCtrl, handles.CompMethod);
% set correction freq range (in Hz)
handles.CorrFrange = [10 10000];
% default normalize status
handles.Normalize = 'on';
% default LowCut options
handles.LowCut = 'off';
handles.LowCutFreq = read_ui_str(handles.LowCutFreqCtrl, 'n');
if strcmpi(handles.LowCut, 'off')
	disable_ui(handles.LowCutFreqText);
	disable_ui(handles.LowCutFreqCtrl);
end
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
%% Ctrl Callbacks
%------------------------------------------------------------------------------
%------------------------------------------------------------------------------
%------------------------------------------------------------------------------

%------------------------------------------------------------------------------
% --- Executes on selection change in CompMethodCtrl.
function CompMethodCtrl_Callback(hObject, eventdata, handles)
	handles.CompMethod = read_ui_val(handles.CompMethodCtrl);
%------------------------------------------------------------------------------


%------------------------------------------------------------------------------
% --- Executes on button press in UpdateSignalCtrl.
%------------------------------------------------------------------------------
function UpdateSignalCtrl_Callback(hObject, eventdata, handles)
	% get mode
	
	switch handles.SignalMode
		
		case 'WAV'
			% load wav file
			[raw, wav.Fs, wav.nbits, wav.opts] = wavread(handles.wavfile);
			
		case 'SYNTH'
			updateSynthFromGui(hObject, handles);
			guidata(hObject, handles);
			synth = handles.synth;
			
			switch synth.type
				case 'tone'
					raw = synmonosine(synth.dur, handles.S.Fs, synth.freq, synth.amp, 0);

				case 'noise'
					raw = synmononoise_fft(synth.dur, handles.S.Fs, synth.fmin, synth.fmax, synth.amp, 0);
					% kludge to scale amplitude properly
					raw = synth.amp * normalize(raw);
					
				case 'sweep'
					raw = synsweep(synth.dur, handles.S.Fs, synth.fmin, synth.fmax, synth.amp, 0);
					
			end
	end
	% apply ramp
	raw = sin2array(raw, synth.ramp, handles.S.Fs);
	% store in handles;
	handles.raw = raw;
	% take fft of raw data
	[handles.fraw, handles.magraw, handles.phiraw] = daqdbfullfft(raw, handles.S.Fs, length(raw));
	guidata(hObject, handles);

	% apply compensation
	switch handles.CompMethod
		case 1
			method = 'ATTEN';
		case 2
			method = 'BOOST'
		case 3
			method = 'COMPRESS'
	end
	if strcmp(handles.LowCut, 'off')
		lowcut = 'off';
	else
		lowcut = handles.LowCutFreq;
	end
	adj = compensate_signal(	raw, ...
										handles.cal.freq, ...
										handles.cal.mag(1, :), ...
										handles.S.Fs, ...
										handles.CorrFrange, ...
										'Method', method, ...
										'Normalize', handles.Normalize, ...
										'Lowcut', lowcut	);
	% store in handles;
	handles.adj = adj;
	% take fft of adj data
	[handles.fadj, handles.magadj, handles.phiadj] = daqdbfullfft(adj, handles.S.Fs, length(adj));
	guidata(hObject, handles);

	updatePlots(hObject, handles);
	guidata(hObject, handles);
%------------------------------------------------------------------------------

%------------------------------------------------------------------------------
function WavSignalButton_Callback(hObject, eventdata, handles)
	newVal = read_ui_val(hObject);
	if newVal == 1
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
		end
% 	else
% 		update_ui_val(handles.WavSignalButton, 1);
	end
	guidata(hObject, handles);
%------------------------------------------------------------------------------

%------------------------------------------------------------------------------
function SynthSignalButton_Callback(hObject, eventdata, handles)
	newVal = read_ui_val(hObject);
	if newVal == 1
		handles.SignalMode = 'SYNTH';
		update_ui_val(handles.WavSignalButton, 0);
	else
		update_ui_val(handles.SynthSignalButton, 1);
	end
	guidata(hObject, handles);
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
	handles.TargetSPL = newVal;
	guidata(hObject, handles);
%------------------------------------------------------------------------------

%------------------------------------------------------------------------------
function PlaySignalCtrl_Callback(hObject, eventdata, handles)
	if (handles.S.Fs == 44100)  && ~isempty(handles.raw) && ~isempty(handles.adj)
		disp('RAW:')
		sound(handles.raw, handles.S.Fs);
		pause(2);
		disp('ADJ:')
		sound(handles.adj, handles.S.Fs);
	end
%------------------------------------------------------------------------------

%------------------------------------------------------------------------------
function NormalizeCtrl_Callback(hObject, eventdata, handles)
	newVal = read_ui_val(handles.NormalizeCtrl);
	if newVal
		handles.Normalize = 'on';
	else
		handles.Normalize = 'off';
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
	
	guidata(hObject, handles);
%------------------------------------------------------------------------------

%-------------------------------------------------------------------------
function loadWavFile(hObject, eventdata, handles)
	[wavfile, wavpath] = uigetfile( '*.wav', ...
												'Load wav file...');
	if wavfile ~=0
		wavdata.datafile = fullfile(wavpath, wavfile);
		[wavdata.raw, wavdata.fs, wavdata.nbits, wavdata.opts] = wavread(wavdata.datafile);
		handles.wavdata = wavdata;
		update_ui_str(handles.FilenameCtrl, wavdata.datafile);
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
function FlatWavMenu_Callback(hObject, eventdata, handles)
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
		ylim([0 100]);
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

%-------------------------------------------------------------------------
%-------------------------------------------------------------------------
function CalMenu_Callback(hObject, eventdata, handles)
%-------------------------------------------------------------------------
%-------------------------------------------------------------------------


%-------------------------------------------------------------------------
%-------------------------------------------------------------------------
function WavMenu_Callback(hObject, eventdata, handles)
%-------------------------------------------------------------------------
%-------------------------------------------------------------------------


%-------------------------------------------------------------------------
%-------------------------------------------------------------------------
function LoadWavMenuItem_Callback(hObject, eventdata, handles)
	loadWavFile(hObject, eventdata, handles);
%-------------------------------------------------------------------------
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
%------------------------------------------------------------------------------
%------------------------------------------------------------------------------
%------------------------------------------------------------------------------

function varargout = temp(varargin)
% TEMP MATLAB code for temp.fig
%      TEMP, by itself, creates a new TEMP or raises the existing
%      singleton*.
%
%      H = TEMP returns the handle to a new TEMP or the handle to
%      the existing singleton*.
%
%      TEMP('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in TEMP.M with the given input arguments.
%
%      TEMP('Property','Value',...) creates a new TEMP or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before temp_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to temp_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help temp

% Last Modified by GUIDE v2.5 11-Oct-2012 18:53:45

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @temp_OpeningFcn, ...
                   'gui_OutputFcn',  @temp_OutputFcn, ...
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
% --- Executes just before temp is made visible.
%------------------------------------------------------------------------------
function temp_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to temp (see VARARGIN)

% Choose default command line output for temp
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);


% define some things
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

guidata(hObject, handles);
updateGuiFromSynth(hObject, handles)
guidata(hObject, handles);
updateSynthFromGui(hObject, handles);
guidata(hObject, handles);

% set compensation method
handles.CompMethod = 1;
% update_ui_str(handles.CompMethodCtrl, {'none'; 'atten'; 'boost'; 'compress'});
update_ui_val(handles.CompMethodCtrl, handles.CompMethod);
guidata(hObject, handles);
update_

% This sets up the initial plot - only do when we are invisible
% so window can get raised using temp.
% if strcmp(get(hObject,'Visible'),'off')
%     plot();
% end

% UIWAIT makes temp wait for user response (see UIRESUME)
% uiwait(handles.figure1);
%------------------------------------------------------------------------------

%------------------------------------------------------------------------------
% --- Outputs from this function are returned to the command line.
%------------------------------------------------------------------------------
function varargout = temp_OutputFcn(hObject, eventdata, handles)
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
% --- Executes on button press in UpdateCtrl.
%------------------------------------------------------------------------------
function UpdateCtrl_Callback(hObject, eventdata, handles)

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

	max(raw)
	handles.raw = raw;
	[handles.fraw, handles.magraw, handles.phiraw] = daqdbfullfft(handles.raw, handles.S.Fs, length(handles.raw));
	guidata(hObject, handles);
	updatePlots(hObject, handles);
	guidata(hObject, handles);
%------------------------------------------------------------------------------

%------------------------------------------------------------------------------
function WavSignalButton_Callback(hObject, eventdata, handles)
	newVal = read_ui_val(hObject);
	if newVal == 1
		handles.SignalMode = 'WAV';
		update_ui_val(handles.SynthSignalButton, 0);
	else
		update_ui_val(handles.WavSignalButton, 1);
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
	title('Raw Signal')
	ylabel('V')
	guidata(hObject, handles);
	
	axes(handles.RawMagAxes)
	plot(0.001*handles.fraw, handles.magraw);
	title('Raw Signal Magnitude')
	xlabel('freq (kHz)')
	ylabel('dB')
	ylim(dblim);
	xlim(freqlim);

	guidata(hObject, handles);
%------------------------------------------------------------------------------

%------------------------------------------------------------------------------
%------------------------------------------------------------------------------
%------------------------------------------------------------------------------
%% MENU Callbacks
%------------------------------------------------------------------------------
%------------------------------------------------------------------------------
%------------------------------------------------------------------------------

%------------------------------------------------------------------------------
function FileMenu_Callback(hObject, eventdata, handles)
%------------------------------------------------------------------------------

%------------------------------------------------------------------------------
function OpenMenuItem_Callback(hObject, eventdata, handles)
	file = uigetfile('*.fig');
	if ~isequal(file, 0)
		 open(file);
	end
%------------------------------------------------------------------------------

%------------------------------------------------------------------------------
function PrintMenuItem_Callback(hObject, eventdata, handles)
	printdlg(handles.figure1)
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

%------------------------------------------------------------------------------
%------------------------------------------------------------------------------
%% Executes during object creation, after setting all properties.
%------------------------------------------------------------------------------
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

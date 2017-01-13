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

% Last Modified by GUIDE v2.5 17-Dec-2015 13:32:48

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
		pdir = ['C:\TytoLogy\TytoLogySettings\' getenv('USERNAME')];
		% development tree
		% pdir = ['C:\Users\sshanbhag\Code\Matlab\TytoLogy\TytoLogySettings\' getenv('USERNAME')];
	else ismac
		pdir = ['~/Work/Code/Matlab/dev/TytoLogy/TytoLogySettings/' ...
							getenv('USER')];
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
	% set output signal to wav (vs. SYNTH) and update GUI
	handles.SignalMode = 'WAV';
	update_ui_val(handles.SynthSignalButton, 0);
	update_ui_val(handles.WavSignalButton, 1);
	% initialize S synth parameter structure
	handles.S = FlatWav_buildS;
	guidata(hObject, handles);
	% initialize tone, noise, sweep structs
	typenum = 1;
	handles.Tone.type = handles.S.Types{typenum};
	for n = 1:handles.S.Nparam(typenum);
		handles.Tone.(handles.S.Param{typenum}{n}) = ...
															handles.S.DefaultVals{typenum}(n);
	end
	typenum = 2;
	handles.Noise.type = handles.S.Types{typenum};
	for n = 1:handles.S.Nparam(typenum);
		handles.Noise.(handles.S.Param{typenum}{n}) = ...
															handles.S.DefaultVals{typenum}(n);
	end
	typenum = 3;
	handles.Sweep.type = handles.S.Types{typenum};
	for n = 1:handles.S.Nparam(typenum);
		handles.Sweep.(handles.S.Param{typenum}{n}) = ...
															handles.S.DefaultVals{typenum}(n);
	end
	% set current synth object to Noise
	handles.synth = handles.Noise;
	handles.SynthIndex = 2;
	handles.SynthType = handles.S.Types{handles.SynthIndex};

	guidata(hObject, handles);
	
	%--------------------------------------------------
	%--------------------------------------------------
	% initialize signals/storage
	%--------------------------------------------------
	%--------------------------------------------------
	handles.raw = [];
	handles.adj = [];
	handles.resp = [];
	handles.rawresp = [];
	handles.rawRMS = [];
	handles.rawdBSPL = [];
	handles.rawfresp = [];
	handles.rawmag = [];
	handles.rawphi = [];
	handles.adjresp = [];
	handles.adjRMS = [];
	handles.adjdBSPL = [];
	handles.adjfresp = [];
	handles.adjmag = [];
	handles.adjphi = [];
	handles.respFs = [];
	guidata(hObject, handles);
	
	%--------------------------------------------------
	%--------------------------------------------------
	% update GUI and synth
	%--------------------------------------------------
	%--------------------------------------------------
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
	set(handles.CompMethodCtrl, 'string', ...
												'none|normalize|atten|boost|compress');
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
    	update_ui_val(handles.LowCutCtrl,0);
    else
		enable_ui(handles.LowCutFreqText);
		enable_ui(handles.LowCutFreqCtrl);
    	update_ui_val(handles.LowCutCtrl,1);
	end
	% target SPL
	handles.TargetSPL = 65;
	update_ui_str(handles.TargetSPLCtrl, handles.TargetSPL);
	% pre-filter range
	handles.PreFilter = 'off';
	handles.PreFilterRange = handles.CorrFrange;
	if strcmpi(handles.PreFilter, 'off')
		update_ui_val(handles.PreFilterCtrl, 1);
	else
		update_ui_val(handles.PreFilterCtrl, 0);
	end
	update_ui_str(handles.PreFilterRangeCtrl, ...
											['[' num2str(handles.PreFilterRange) ']']);
	disable_ui(handles.PreFilterRangeCtrl);
	guidata(hObject, handles);
	% post-filter range
	handles.PostFilter = 'off';
	handles.PostFilterRange = handles.CorrFrange;
	if strcmpi(handles.PostFilter, 'off')
		update_ui_val(handles.PostFilterCtrl, 1);
	else
		update_ui_val(handles.PostFilterCtrl, 0);
	end
	update_ui_str(handles.PostFilterRangeCtrl, ...
											['[' num2str(handles.PostFilterRange) ']']);
	if strcmpi(handles.PostFilter, 'off')
		disable_ui(handles.PostFilterRangeCtrl);
	else
		enable_ui(handles.PostFilterRangeCtrl);
	end

	guidata(hObject, handles);
	% range limit
	handles.RangeLimit = 'off';
	if strcmpi(handles.RangeLimit, 'off')
		update_ui_val(handles.RangeLimitCtrl, 1);
	else
		update_ui_val(handles.RangeLimitCtrl, 0);
	end
	guidata(hObject, handles);
	% CorrectionLimit
	handles.CorrectionLimit = 'off';
	handles.CorrectionLimitValue = 5;
	if strcmpi(handles.CorrectionLimit, 'off')
		update_ui_val(handles.CorrectionLimitCtrl, 1);
	else
		update_ui_val(handles.CorrectionLimitCtrl, 0);
	end
	update_ui_str(handles.CorrectionLimitValCtrl, handles.CorrectionLimitValue);
	if strcmpi(handles.CorrectionLimit, 'off')
		disable_ui(handles.CorrectionLimitValCtrl);
		disable_ui(handles.CorrectionLimitText);
	else
		enable_ui(handles.CorrectionLimitValCtrl);
		enable_ui(handles.CorrectionLimitText);
	end
	guidata(hObject, handles);
	% SmoothEdges
	handles.SmoothEdges= 'off';
	handles.SmoothEdgesValue = [1 5];
	if strcmpi(handles.SmoothEdges, 'off')
		update_ui_val(handles.SmoothEdgesCtrl, 1);
	else
		update_ui_val(handles.SmoothEdgesCtrl, 0);
	end
	update_ui_str(handles.SmoothEdgesValCtrl, ...
											['[' num2str(handles.SmoothEdgesValue) ']']);
	if strcmpi(handles.SmoothEdges, 'off')
		disable_ui(handles.SmoothEdgesValCtrl);
	else
		enable_ui(handles.SmoothEdgesValCtrl);
	end
	guidata(hObject, handles);

	%--------------------------------------------------
	%--------------------------------------------------
	% spectrum settings
	%--------------------------------------------------
	%--------------------------------------------------
	handles.PlotSpectrum = 'off';
	set(handles.PlotSpectrumMenuItem, 'Checked', handles.PlotSpectrum);
	handles.SpectrumWindow = 512;
	handles.ColorMap = 'gray';
	guidata(hObject, handles);
	
	
	%--------------------------------------------------
	%--------------------------------------------------
	% input filter settings
	%--------------------------------------------------
	handles.HPFc = 100; 
	handles.LPFc = 100000;
	handles.FilterOrder = 3;
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
	plot(	handles.CalibrationAxes, ...
			0.001*handles.cal.freq, ...
			handles.cal.mag(1, :), '.-');
	ylim(handles.CalibrationAxes, [0 100]);
	
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
	handles.OutputDevice = 'NIDAQ';
	handles.IOpause = 0.5;
	update_ui_val(handles.SoundCardButton, 0);
	update_ui_val(handles.NIDAQButton, 1);
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
		
	%--------------------------------------------------
	%--------------------------------------------------
	% analysis settings
	%--------------------------------------------------
	%--------------------------------------------------
	% set Analysis window
	handles.Awindow = [0 handles.S.DefaultVals{typenum}(1)];
	update_ui_str(handles.AnalysisStartCtrl, handles.Awindow(1));
	update_ui_str(handles.AnalysisEndCtrl, handles.Awindow(2));
	% window for finding RMS peak, in milliseconds
	handles.PeakRMSWindow = 5;
	handles.dBPlot = 0;
	handles.dBFigure = [];
	guidata(hObject, handles);
%------------------------------------------------------------------------------
%------------------------------------------------------------------------------
%******************************************************************************
%******************************************************************************
%******************************************************************************


%******************************************************************************
%******************************************************************************
%******************************************************************************
%------------------------------------------------------------------------------
% --- Outputs from this function are returned to the command line.
%------------------------------------------------------------------------------
function varargout = FlatWav_OutputFcn(hObject, eventdata, handles) %#ok<*INUSL>
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
% ACTION BUTTON CONTROL Callbacks (on upper left side of GUI)
%******************************************************************************
%******************************************************************************
%******************************************************************************

%------------------------------------------------------------------------------
% --- Executes on button press in UpdateSignalCtrl.
%------------------------------------------------------------------------------
function UpdateSignalCtrl_Callback(hObject, eventdata, handles) %#ok<*DEFNU>
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
			if verLessThan('matlab', 'R2015b')
				[wavdata.raw, wavdata.Fs, wavdata.nbits, wavdata.opts] = ...
									wavread(handles.wavdata.datafile); %#ok<DWVRD>

			else
				tmp = audioinfo(handles.wavdata.datafile);
				wavdata.Fs = tmp.SampleRate;
				wavdata.nbits = tmp.BitsPerSample';
				wavdata.opts = tmp;
				wavdata.raw = audioread(handles.wavdata.datafile);
			end
			wavdata.datafile = handles.wavdata.datafile;
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
	% next, check the calibration curve and smooth if necessary
	%--------------------------------------------------------------------
	if read_ui_val(handles.SmoothCalCtrl)
		% get value of smooth method
		smoothmethod = read_ui_val(handles.CalSmoothMethodCtrl);
		switch(smoothmethod)
			% moving-window average
			case 1
				% smooth the calibration mag data, store in mag_smooth
				window_len = read_ui_str(handles.SmoothVal1Ctrl, 'n');
				handles.mag_smooth = smooth_calibration_data(...
																		smoothmethod, ...
																		handles.cal, ...
																		window_len);
			% Savitzky-Golay filter
			case 2
				% smooth the calibration mag data, store in mag_smooth
				sg_order = read_ui_str(handles.SmoothVal1Ctrl, 'n');
				sg_framesize = read_ui_str(handles.SmoothVal2Ctrl, 'n');
				handles.mag_smooth = smooth_calibration_data( ...	
																		smoothmethod, ...
																		handles.cal, ...
																		sg_order, ...
																		sg_framesize);
		end
		guidata(hObject, handles);
		mags = handles.mag_smooth;
	else
		% use unsmoothed mags
		mags = handles.cal.mag;
	end
	
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
	if strcmpi(handles.LowCut, 'off')
		lowcut = 'off';
	else
		lowcut = handles.LowCutFreq;
	end
	% check pre-filter
	if strcmpi(handles.PreFilter, 'off')
		prefilter = 'off';
	else
		prefilter = handles.PreFilterRange;
	end
	% check post-filter
	if strcmpi(handles.PostFilter, 'off')
		postfilter = 'off';
	else
		postfilter = handles.PostFilterRange;
	end
	% check rangelimit
	if strcmpi(handles.RangeLimit, 'off')
		rangelimit = 'off';
	else
		rangelimit = 'on';
	end
	% check Correction Limit
	if strcmpi(handles.CorrectionLimit, 'off')
		corrlimit = 'off';
	else
		corrlimit = handles.CorrectionLimitValue;
	end
	% check SmoothEdges
	if strcmpi(handles.SmoothEdges, 'off')
		smoothedges = 'off';
	else
		smoothedges = handles.SmoothEdgesValue;
	end
	
	% compensate signal by method
	if strcmpi(method, 'NONE')
		% no compensation, just copy raw to adj
		handles.adj = handles.raw;
	
	elseif strcmpi(method, 'NORMALIZE')
		% normalize only, using Normalize Value
		handles.adj = handles.NormalizeValue * normalize(handles.raw);
		
	elseif strcmpi(handles.Normalize, 'off')
		% kludge to deal with Normalize OFF vs. on
		[handles.adj, tmp, handles.compcurve] = compensate_signal(	...
											handles.raw, ...
											handles.cal.freq, ...
											mags(1, :), ...
											handles.S.Fs, ...
											handles.CorrFrange, ...
											'Method', method, ...
											'Normalize', 'off', ...
											'Lowcut', lowcut, ...
											'Level', handles.TargetSPL, ...
											'Prefilter', prefilter, ...
											'Postfilter', postfilter, ...
											'Rangelimit', rangelimit, ...
											'Corrlimit', corrlimit, ...
											'SmoothEdges', smoothedges	); %#ok<*ASGLU>
		clear tmp

	else
		% use normalization
		[handles.adj, tmp, handles.compcurve] = compensate_signal(	...
											handles.raw, ...
											handles.cal.freq, ...
											mags(1, :), ...
											handles.S.Fs, ...
											handles.CorrFrange, ...
											'Method', method, ...
											'Normalize', handles.NormalizeValue, ...
											'Lowcut', lowcut, ...
											'Level', handles.TargetSPL, ...
											'Prefilter', prefilter, ...
											'Postfilter', postfilter, ...
											'Rangelimit', rangelimit, ...
											'Corrlimit', corrlimit, ...
											'SmoothEdges', smoothedges	);
		clear tmp
	end
	
	if strcmpi(handles.SignalMode, 'SYNTH')
		% apply ramp to stimuli if it was synthesized
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
	updatePlots(hObject, eventdata, handles);
	% save variables in workspace
	guidata(hObject, handles);
	disp('....updating complete');
%------------------------------------------------------------------------------

%------------------------------------------------------------------------------
function PlayRawSignalCtrl_Callback(hObject, eventdata, handles)
%------------------------------------------------------------------------------
	[hObject, handles] = PlaySignal(hObject, eventdata, handles);
	guidata(hObject, handles);
%------------------------------------------------------------------------------

%------------------------------------------------------------------------------
function PlayAdjSignalCtrl_Callback(hObject, eventdata, handles)
%------------------------------------------------------------------------------
	[hObject, handles] = PlaySignal(hObject, eventdata, handles);
	guidata(hObject, handles);
%------------------------------------------------------------------------------

%------------------------------------------------------------------------------
function SaveSoundCtrl_Callback(hObject, eventdata, handles)
	% use the menu item callback
	SaveAdjSignalMenuItem_Callback(hObject, eventdata, handles);
%------------------------------------------------------------------------------


%******************************************************************************
%******************************************************************************
%******************************************************************************

%{
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
		ylim(handles.CalibrationAxes,...
					[	(0.9 * min(handles.cal.mag(1, :))) ...
					(1.1 * max(handles.cal.mag(1, :)))		]);
		grid(handles.CalibrationAxes, 'on');
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
			% smooth the calibration mag data, store in mag_smooth
			window_len = read_ui_str(handles.SmoothVal1Ctrl, 'n');
			mag_smooth = smooth_calibration_data(	smoothmethod, ...
																handles.cal, ...
																window_len);
		case 2
			% Savitzky-Golay filter
			enable_ui(handles.SmoothVal1Text);
			update_ui_str(handles.SmoothVal1Text, 'Order');
			enable_ui(handles.SmoothVal1Ctrl);
			enable_ui(handles.SmoothVal2Text);
			update_ui_str(handles.SmoothVal2Text, 'Frame Size');
			enable_ui(handles.SmoothVal2Ctrl);
			% smooth the calibration mag data, store in mag_smooth
			sg_order = read_ui_str(handles.SmoothVal1Ctrl, 'n');
			sg_framesize = read_ui_str(handles.SmoothVal2Ctrl, 'n');
			mag_smooth = smooth_calibration_data(	smoothmethod, ...
																handles.cal, ...
																sg_order, ...
																sg_framesize);
	end
	
	% plot the smoothed data
	plot(handles.CalibrationAxes, ...
				0.001*handles.cal.freq, ...
				mag_smooth(1, :), 'b.-');
	ylim(handles.CalibrationAxes, ...
		[(0.9 * min(mag_smooth(1, :))) (1.1 * max(mag_smooth(1, :)))]);
	grid(handles.CalibrationAxes, 'on');
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
%}

%------------------------------------------------------------------------------
function LoadCalCtrl_Callback(hObject, eventdata, handles)
	% use the menu item callback
	LoadCalMenuItem_Callback(hObject, eventdata, handles)
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
	param = handles.S.Param{handles.SynthIndex}{tagnum}; %#ok<FNDSB>
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
		% load wav file
		if verLessThan('matlab', 'R2015b')
			[wavdata.raw, wavdata.Fs, wavdata.nbits, wavdata.opts] = ...
													wavread(wavdata.datafile); %#ok<DWVRD>
		else
			tmp = audioinfo(wavdata.datafile);
			wavdata.Fs = tmp.SampleRate;
			wavdata.nbits = tmp.BitsPerSample';
			wavdata.opts = tmp;
			wavdata.raw = audioread(wavdata.datafile);
		end
		handles.wavdata = wavdata;
		guidata(hObject, handles);
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
%------------------------------------------------------------------------------

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
function CompMethodCtrl_Callback(hObject, eventdata, handles)
%-------------------------------------------------------------------
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
		NormalizePeakCtrl_Callback(hObject, eventdata, handles);
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
%-------------------------------------------------------------------
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
%-------------------------------------------------------------------
	handles.LowCutFreq = read_ui_str(hObject, 'n');
	guidata(hObject, handles);
%------------------------------------------------------------------------------

%------------------------------------------------------------------------------
function CorrFminCtrl_Callback(hObject, eventdata, handles)
%-------------------------------------------------------------------
	newVal = read_ui_str(hObject, 'n');
	% NEED CHECKS!
	handles.CorrFrange(1) = newVal;
	guidata(hObject, handles);
%------------------------------------------------------------------------------

%------------------------------------------------------------------------------
function CorrFmaxCtrl_Callback(hObject, eventdata, handles)
%-------------------------------------------------------------------
	newVal = read_ui_str(hObject, 'n');
	% NEED CHECKS!
	handles.CorrFrange(2) = newVal;
	guidata(hObject, handles);
%------------------------------------------------------------------------------

%------------------------------------------------------------------------------
function PreFilterCtrl_Callback(hObject, eventdata, handles)
%-------------------------------------------------------------------
	newVal = read_ui_val(handles.PreFilterCtrl);
	if newVal
		handles.PreFilter = 'on';
		enable_ui(handles.PreFilterRangeCtrl);
		PreFilterRangeCtrl_Callback(hObject, eventdata, handles);
	else
		handles.PreFilter = 'off';
		disable_ui(handles.PreFilterRangeCtrl);
	end
	guidata(hObject, handles)
%------------------------------------------------------------------------------

%------------------------------------------------------------------------------
function PreFilterRangeCtrl_Callback(hObject, eventdata, handles)
%-------------------------------------------------------------------
	try
		tmpstr = read_ui_str(handles.PreFilterRangeCtrl);
		newVal = eval(tmpstr);
		if newVal(1) > newVal(2)
			errordlg('Bad filter range', 'Flat Wav Error');
			update_ui_str(hObject, ['[' num2str(handles.PreFilterRange) ']']);
		else
			handles.PreFilterRange = newVal;
			guidata(hObject, handles);
		end
	catch err %#ok<NASGU>
		fprintf('error in PreFilterRangeCtrl_Callback\n');
		disp err
	end
%------------------------------------------------------------------------------

%------------------------------------------------------------------------------
function PostFilterCtrl_Callback(hObject, eventdata, handles)
	newVal = read_ui_val(handles.PostFilterCtrl);
	if newVal
		handles.PostFilter = 'on';
		enable_ui(handles.PostFilterRangeCtrl);
		PostFilterRangeCtrl_Callback(hObject, eventdata, handles);
	else
		handles.PostFilter = 'off';
		disable_ui(handles.PostFilterRangeCtrl);
	end
	guidata(hObject, handles)
%------------------------------------------------------------------------------

%------------------------------------------------------------------------------
function PostFilterRangeCtrl_Callback(hObject, eventdata, handles)
%-------------------------------------------------------------------
	try
		tmpstr = read_ui_str(handles.PostFilterRangeCtrl);
		newVal = eval(tmpstr);
		if newVal(1) > newVal(2)
			errordlg('Bad filter range', 'Flat Wav Error');
			update_ui_str(hObject, ['[' num2str(handles.PostFilterRange) ']']);
		else
			handles.PostFilterRange = newVal;
			guidata(hObject, handles);
		end
	catch err %#ok<NASGU>
		fprintf('error in PostFilterRangeCtrl_Callback\n');
		disp err
	end
%------------------------------------------------------------------------------

%------------------------------------------------------------------------------
function RangeLimitCtrl_Callback(hObject, eventdata, handles)
%-------------------------------------------------------------------
	newVal = read_ui_val(hObject);
	if newVal
		handles.RangeLimit = 'on';
	else
		handles.RangeLimit = 'off';
	end
	guidata(hObject, handles);
%------------------------------------------------------------------------------

%------------------------------------------------------------------------------
function CorrectionLimitCtrl_Callback(hObject, eventdata, handles)
%-------------------------------------------------------------------
	newVal = read_ui_val(hObject);
	if newVal
		handles.CorrectionLimit = 'on';
		enable_ui(handles.CorrectionLimitValCtrl);
	else
		handles.CorrectionLimit = 'off';
		disable_ui(handles.CorrectionLimitValCtrl);
	end
	guidata(hObject, handles)
%------------------------------------------------------------------------------

%------------------------------------------------------------------------------
function CorrectionLimitValCtrl_Callback(hObject, eventdata, handles)
%-------------------------------------------------------------------
	try
		newVal = read_ui_str(hObject, 'n');
		if newVal <= 0
			errordlg('Bad Correction Limit', 'Flat Wav Error');
			update_ui_str(hObject, handles.CorrectionLimitValue);
		else
			handles.CorrectionLimitValue = newVal;
			guidata(hObject, handles);
		end
	catch err %#ok<NASGU>
		fprintf('error in CorrectionLimitValCtrl_Callback\n');
		disp err
	end
%------------------------------------------------------------------------------

%------------------------------------------------------------------------------
function SmoothEdgesCtrl_Callback(hObject, eventdata, handles)
%-------------------------------------------------------------------
	newVal = read_ui_val(hObject);
	if newVal
		handles.SmoothEdges = 'on';
		enable_ui(handles.SmoothEdgesValCtrl);
	else
		handles.SmoothEdges = 'off';
		disable_ui(handles.SmoothEdgesValCtrl);
	end
	guidata(hObject, handles)
%------------------------------------------------------------------------------

%------------------------------------------------------------------------------
function SmoothEdgesValCtrl_Callback(hObject, eventdata, handles)
%-------------------------------------------------------------------
	try
		tmpstr = read_ui_str(hObject);
		newVal = eval(tmpstr);
		if newVal(1) > newVal(2)
			errordlg('Bad SmoothEdges values', 'Flat Wav Error');
			update_ui_str(hObject, ['[' num2str(handles.SmoothEdgesValue) ']']);
		else
			handles.SmoothEdgesValue = newVal;
			guidata(hObject, handles);
		end
	catch err %#ok<NASGU>
		fprintf('error in SmoothEdgesValCtrl_Callback\n');
		disp err
	end
%------------------------------------------------------------------------------

%------------------------------------------------------------------------------
function PlotCompensationCurveCtrl_Callback(hObject, eventdata, handles)
	if ~isfield(handles, 'cal')
		return
	elseif isempty(handles.cal)
		return
	elseif ~isfield(handles, 'compcurve')
		return
	elseif isempty(handles.compcurve)
		return
	end
	figure
	plot(0.001 * handles.cal.freq, handles.compcurve, 'k.-');
	xlabel('Frequency (kHz)')
	ylabel('dB');
%------------------------------------------------------------------------------
%------------------------------------------------------------------------------
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
function update_analysis(hObject, eventdata, handles)
	if ~isempty(handles.respFs)
		rmsbins = ms2samples(handles.PeakRMSWindow, handles.respFs); %#ok<NASGU>
	else
		fprintf('\n\nupdate_analysis: handles.respFs is empty!\n\n');
		return
	end
	%--------------------------------------------------
	% update analysis window
	%--------------------------------------------------
	% find max value of raw response
	if ~isempty(handles.rawresp) 
		% check if analysis window is beyond length of signal
		if  ms2samples(handles.Awindow(2), handles.respFs) > ...
										length(handles.rawresp)
			% if so, reset to duration of signal
			handles.Awindow(2) = floor(bin2ms(length(handles.rawresp), ...
													handles.respFs));
			update_ui_str(handles.AnalysisEndCtrl, handles.Awindow(2));
			guidata(hObject, handles);
			fprintf('warning: Analysis End > length of raw signal!!!!');
		end
		% find bins for analysis
		bin = ms2samples(handles.Awindow, handles.respFs);
		if bin(1) == 0
			bin(1) = 1;
		end
		%--------------------------------------------------
		% analyze data
		%--------------------------------------------------
		% compute RMS
		resp_RMS = rms(handles.rawresp(bin(1):bin(2)));
		% compute dB SPL
		 handles.rawdBSPL = dbspl(handles.VtoPa*resp_RMS);
		% update display
		dbtext = sprintf('Raw dB SPL: %.2f  [%d - %d]\n', ...
															handles.rawdBSPL, ...
															handles.Awindow(1), ...
															handles.Awindow(2));
		fprintf('%s\n', dbtext);
		update_ui_str(handles.RawdBText, dbtext);
		show_uictrl(handles.RawdBText);
	end
	%--------------------------------------------------
	% update analysis window
	%--------------------------------------------------
	% find max value of raw response
	if ~isempty(handles.adjresp)
		% check if analysis window is beyond length of signal
		if  ms2samples(handles.Awindow(2), handles.respFs) > ...
							length(handles.adjresp)
			% if so, reset to duration of signal
			handles.Awindow(2) = floor(bin2ms(length(handles.adjresp), ...
								handles.respFs));
			update_ui_str(handles.AnalysisEndCtrl, handles.Awindow(2));
			guidata(hObject, handles);
			fprintf('warning: Analysis End > length of adj signal!!!!');
		end
		% find bins for analysis
		bin = ms2samples(handles.Awindow, handles.respFs);
		if bin(1) == 0
			bin(1) = 1;
		end
		%--------------------------------------------------
		% analyze data
		%--------------------------------------------------
		% compute RMS
		resp_RMS = rms(handles.adjresp(bin(1):bin(2)));
		% compute dB SPL
		 handles.adjdBSPL = dbspl(handles.VtoPa*resp_RMS);
		% update display
		dbtext = sprintf('Raw dB SPL: %.2f  [%d - %d]\n', ...
															handles.adjdBSPL, ...
															handles.Awindow(1), ...
															handles.Awindow(2));
		fprintf('%s\n', dbtext);
		update_ui_str(handles.AdjdBText, dbtext);
		show_uictrl(handles.AdjdBText);
	end
	guidata(hObject, handles);
%------------------------------------------------------------------------------

%------------------------------------------------------------------------------
function AnalysisStartCtrl_Callback(hObject, eventdata, handles)
%-------------------------------------------------------------------
	newVal = read_ui_str(hObject, 'n');
	if ~isnumeric(newVal)
		update_ui_str(hObject, handles.handles.Awindow(1));
		errordlg('AnalysisStart: Value must be a number', 'Flat Wav Error');
	elseif newVal <= 0
		update_ui_str(hObject, handles.handles.Awindow(1));
		errordlg('AnalysisStart: Value must > 0', 'Flat Wav Error');
	end
	handles.Awindow(1) = newVal;	
	guidata(hObject, handles);
	update_analysis(hObject, eventdata, handles);
%------------------------------------------------------------------------------

%------------------------------------------------------------------------------
function AnalysisEndCtrl_Callback(hObject, eventdata, handles)
%-------------------------------------------------------------------
	newVal = read_ui_str(hObject, 'n');
	if ~isnumeric(newVal)
		update_ui_str(hObject, handles.handles.Awindow(2));
		errordlg('AnalysisEnd: Value must be a number', 'Flat Wav Error');
	elseif newVal <= 0
		update_ui_str(hObject, handles.handles.Awindow(2));
		errordlg('AnalysisEnd: Value must > 0', 'Flat Wav Error');
	else
		handles.Awindow(2) = newVal;
	end
	guidata(hObject, handles);
	update_analysis(hObject, eventdata, handles)
%------------------------------------------------------------------------------

%------------------------------------------------------------------------------
function FindPeakCtrl_Callback(hObject, eventdata, handles)
%-------------------------------------------------------------------
	if ~isempty(handles.respFs)
		rmsbins = ms2samples(handles.PeakRMSWindow, handles.respFs);
	else
		fprintf('\n\nFindPeakCtrl_Callback: handles.respFs is empty!\n\n');
		return
	end
	% find max value of raw response
	if ~isempty(handles.rawresp)
		tmprms = block_rms(handles.rawresp, rmsbins);
		% find peak and peak index of rms values
		[handles.rawrespmax.val, handles.rawrespmax.indx] = max(tmprms);
		clear temprms
		% compute peak dB SPL
		rawdBSPL = dbspl(handles.VtoPa*handles.rawrespmax.val);
		% update display
		dbtext = sprintf('Raw dB SPL: %.2f  [%d - %d]\n', ...
															rawdBSPL, ...
															handles.Awindow(1), ...
															handles.Awindow(2));
		fprintf('%s\n', dbtext);
		update_ui_str(handles.RawdBText, dbtext);
		show_uictrl(handles.RawdBText);
		% find max point (in milliseconds)
		xval = rmsbins * handles.rawrespmax.indx + (rmsbins ./ 2);
		xval = fix(bin2ms(xval, handles.respFs));
		update_ui_str(handles.PeakTimeRawCtrl, xval);
		guidata(hObject, handles);
	end
	% find max value of adj response
	if ~isempty(handles.adjresp)
		tmprms = block_rms(handles.adjresp, rmsbins);
		% find peak and peak index of rms values
		[handles.adjrespmax.val, handles.adjrespmax.indx] = max(tmprms);
		clear temprms
		% compute peak dB SPL
		adjdBSPL = dbspl(handles.VtoPa*handles.rawrespmax.val);
		% update display
		dbtext = sprintf('Raw dB SPL: %.2f  [%d - %d]\n', ...
															adjdBSPL, ...
															handles.Awindow(1), ...
															handles.Awindow(2));
		fprintf('%s\n', dbtext);
		update_ui_str(handles.AdjdBText, dbtext);
		show_uictrl(handles.AdjdBText);
		% find max point (in milliseconds)
		xval = rmsbins * handles.adjrespmax.indx + (rmsbins ./ 2);
		xval = fix(bin2ms(xval, handles.respFs));
		update_ui_str(handles.PeakTimeAdjCtrl, xval);
		guidata(hObject, handles);
	end
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
		update_ui_str(handles.(handles.S.TextTags{n}), ...
															[handles.S.Text{sindx}{n} ':']);
		update_ui_str(handles.(handles.S.CtrlTags{n}), ...
											handles.synth.(handles.S.Param{sindx}{n}));
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
function updatePlots(hObject, eventdata, handles)
	% plotting limits
	% limits for Mag and Phase plots
	dbmax = 0;
	dbmin = min([min(handles.magraw) min(handles.magadj)]);
	dblim = [dbmin dbmax];
	freqlim = 0.001*[0 handles.S.Fs/2];
	
	% update raw plots
	% plot raw signal
	tvec = 1000 * (0:(length(handles.raw)-1)) ./ handles.S.Fs;
	plot(handles.RawSignalAxes, tvec, handles.raw)
	title(handles.RawSignalAxes, 'Signal (V)')
	ylabel(handles.RawSignalAxes, 'Raw', 'Color', 'b')
	set(handles.RawSignalAxes, 'XTickLabel', []);
	xlim(handles.RawSignalAxes, [min(tvec) max(tvec)])
	% get ticks
	time_ticks = get(handles.RawSignalAxes, 'XTick');
	
	% 	plot raw spectrogram 
	if strcmpi(handles.PlotSpectrum, 'on')
		% plot raw magnitude
		plot(handles.RawMagAxes, 0.001*handles.fraw, handles.magraw);
		title(handles.RawMagAxes, 'Magnitude (dB)')
		ylim(handles.RawMagAxes, dblim);
		xlim(handles.RawMagAxes, freqlim);
		set(handles.RawMagAxes, 'XTickLabel', []);

		% plot raw phase
		plot(handles.RawPhaseAxes, 0.001*handles.fraw, unwrap(handles.phiraw));
		title(handles.RawPhaseAxes, 'Phase (rad)')
		xlim(handles.RawPhaseAxes, freqlim);
		set(handles.RawPhaseAxes, 'XTickLabel', []);

		% plot raw spectrogram
		[S, F, T, P] = spectrogram(	handles.raw, ...
												handles.SpectrumWindow, ...
												[], ...
												handles.SpectrumWindow, ...
												handles.S.Fs	);
		P = 20*log10(P);
		P(P == -Inf) = min(min(P(P ~= -Inf)));	
		surf(handles.RawSpectrumAxes, 1000*T, 0.001*F, P, 'edgecolor', 'none');
		xlim(handles.RawSpectrumAxes, [min(tvec) max(tvec)])
		ylim(handles.RawSpectrumAxes, freqlim);
		set(handles.RawSpectrumAxes, 'XTick', time_ticks)
		view(handles.RawSpectrumAxes, 0, 90);
		title(handles.RawSpectrumAxes, 'Time vs. Freq (kHz) vs. dB')
		set(handles.RawSpectrumAxes, 'XTickLabel', []);
		colormap(handles.RawSpectrumAxes, handles.ColorMap)
	else
		cla(handles.RawMagAxes);
		cla(handles.RawPhaseAxes);
		cla(handles.RawSpectrumAxes);
	end
	guidata(hObject, handles)
	
	% Update adj plots
	% plot adj signal
	tvec = 1000 * (0:(length(handles.adj)-1)) ./ handles.S.Fs;
	plot(handles.AdjSignalAxes, tvec, handles.adj, 'r')
	xlim(handles.AdjSignalAxes, [min(tvec) max(tvec)])
	ylabel(handles.AdjSignalAxes, 'Adj', 'Color', 'r')
	xlabel(handles.AdjSignalAxes, 'time (ms)')
	
	if strcmpi(handles.PlotSpectrum, 'on')
		% plot adj mag spectrum
		plot(handles.AdjMagAxes, 0.001*handles.fadj, handles.magadj, 'r');
		ylim(handles.AdjMagAxes, dblim);
		xlim(handles.AdjMagAxes, freqlim);
		xlabel(handles.AdjMagAxes, 'freq (kHz)');

		% plot adj phase spectrum
		plot(handles.AdjPhaseAxes, 0.001*handles.fadj, unwrap(handles.phiadj), 'r');
		xlim(handles.AdjPhaseAxes, freqlim);
		xlabel(handles.AdjPhaseAxes, 'freq (kHz)');		
		
		% plot adj spectrogram
		[S, F, T, P] = spectrogram(	handles.adj, ...
												handles.SpectrumWindow, ...
												[], ...
												handles.SpectrumWindow, ...
												handles.S.Fs	);
		P = 20*log10(P);
		P(P == -Inf) = min(min(P(P ~= -Inf)));	
		surf(handles.AdjSpectrumAxes, 1000*T, 0.001*F, P, 'edgecolor', 'none');
		xlim(handles.AdjSpectrumAxes, [min(tvec) max(tvec)])
		ylim(handles.AdjSpectrumAxes, freqlim);
		set(handles.AdjSpectrumAxes, 'XTick', time_ticks)	
		view(handles.AdjSpectrumAxes, 0, 90);
		xlabel(handles.AdjSpectrumAxes, 'Time (ms)')
		colormap(handles.AdjSpectrumAxes, handles.ColorMap)
	else
		cla(handles.AdjMagAxes);
		cla(handles.AdjPhaseAxes);
		cla(handles.AdjSpectrumAxes);
	end
	guidata(hObject, handles);
	
% 	dBPlotCtrl_Callback(hObject, eventdata, handles);
% 	guidata(hObject, handles);
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

function varargout = PlaySignal(hObject, eventdata, handles)
%------------------------------------------------------------------------------
% PlaySignal
%------------------------------------------------------------------------------
% sets up NI data acquisition toolbox parameters and plays signal
%------------------------------------------------------------------------------
	%-----------------------------------------------------------------------
	%-----------------------------------------------------------------------
	% Check Inputs
	%-----------------------------------------------------------------------
	%-----------------------------------------------------------------------
	%--------------------------------------------------
	% find who called us
	%--------------------------------------------------
	ButtonID = read_ui_str(hObject);
	% check to make sure output signal isn't crazy
	if strcmpi(ButtonID, 'Play Raw') 
		if isempty(handles.raw)	
			warndlg('RAW signal empty!');
			return
		end
	end
	if strcmpi(ButtonID, 'Play Adj')
		if isempty(handles.adj)
			warndlg('ADJ signal empty!');
			return
		end
	end
	%-----------------------------------------------------------------------
	%-----------------------------------------------------------------------
	% setup output figures
	%-----------------------------------------------------------------------
	%-----------------------------------------------------------------------
	% if figure has not been created, do so
	if isempty(handles.IOfigure) || ~ishandle(handles.IOfigure)
		handles.IOfigure = figure;
		set(handles.IOfigure, 'Position', [10 235 973 500]);
		% create subplot axes
		% 	subplot('Position',[left bottom width height]) creates an axes at the
		% 	position specified by a four-element vector. left, bottom, width, and
		% 	height are in normalized coordinates in the range from 0.0 to 1.0
		pw = 0.19;		% plot width
		ph = 0.375;		% plot height
		x1 = 0.052;
		x2 = 0.3;
		x3 = 0.55;
		x4 = 0.8;
		y1 = 0.55;
		y2 = 0.09;
		handles.P.rsig =	subplot('Position', [x1	y1	pw	ph]);
		handles.P.rmag =	subplot('Position', [x2	y1	pw	ph]);
		handles.P.rphi =	subplot('Position', [x3	y1	pw	ph]);
		handles.P.rspec =	subplot('Position', [x4	y1	pw	ph]);
		handles.P.asig =	subplot('Position', [x1	y2	pw	ph]);
		handles.P.amag =	subplot('Position', [x2	y2	pw	ph]);
		handles.P.aphi =	subplot('Position', [x3	y2	pw	ph]);
		handles.P.aspec =	subplot('Position', [x4	y2	pw	ph]);
		guidata(hObject, handles);
	else
		% otherwise, make it active
		figure(handles.IOfigure);
	end

	%-----------------------------------------------------------------------
	%-----------------------------------------------------------------------
	% play sound out of selected output device
	%-----------------------------------------------------------------------
	%-----------------------------------------------------------------------
	% SOUND CARD OUTPUT
	if strcmpi(handles.OutputDevice, 'WINSOUND')
		if (handles.S.Fs == 44100)
			if strcmpi(ButtonID, 'Play Raw') && ~isempty(handles.raw)
				% play raw sound
				disable_ui(hObject);
				pause(0.5);
				sound(sin2array(handles.raw, 1, handles.S.Fs), handles.S.Fs);
				pause(1);
				enable_ui(hObject);
				update_ui_str(handles.RawdBText, sprintf('Raw dB SPL: ---'));
				hide_uictrl(handles.RawdBText);
				hide_uictrl(handles.AdjdBText);
			elseif strcmpi(ButtonID, 'Play Adj') && ~isempty(handles.adj)
				% play adj sound
				disable_ui(hObject);
				pause(0.5);
				sound(sin2array(handles.adj, 1, handles.S.Fs), handles.S.Fs);
				pause(1);
				enable_ui(hObject);
				update_ui_str(handles.AdjdBText, sprintf('Adj dB SPL: ---'));
				hide_uictrl(handles.RawdBText);
				hide_uictrl(handles.AdjdBText);		
			end
		else
			errordlg(sprintf('Invalid sample rate %.2f for WINSOUND', handles.S.Fs), ...
									'FlatWav Error');
		end
	% NIDAQ OUTPUT
	elseif strcmpi(handles.OutputDevice, 'NIDAQ')
		% play selected sound, get response
		if strcmpi(ButtonID, 'Play Raw') && ~isempty(handles.raw)
			% play raw sound
			disable_ui(hObject);
			[resp, Fs] = NIplaysignal(hObject, handles, handles.raw);
			enable_ui(hObject);		
		elseif strcmpi(ButtonID, 'Play Adj') && ~isempty(handles.adj)
			% play adj sound
			disable_ui(hObject);
			[resp, Fs] = NIplaysignal(hObject, handles, handles.adj);
			show_uictrl(handles.AdjdBText);		
			enable_ui(hObject);
		end
	% ERROR
	else
		errordlg(sprintf('unknown io device %s', handles.OutputDevice), 'FlatWav Error');
	end
	handles.respFs = Fs;
	guidata(hObject, handles);

	%-----------------------------------------------------------------------
	%-----------------------------------------------------------------------
	% process response data (if NIDAQ selected)
	%-----------------------------------------------------------------------
	%-----------------------------------------------------------------------
	if strcmpi(handles.OutputDevice, 'NIDAQ')
		%-----------------------------------------------------------------------
		% update bandpass filter for processing the data
		%-----------------------------------------------------------------------
		% Nyquist frequency
		fnyq = Fs / 2;
		% passband definition
		fband = [handles.HPFc handles.LPFc] ./ fnyq;
		% filter coefficients using a butterworth bandpass filter
		[handles.fcoeffb, handles.fcoeffa] = butter(handles.FilterOrder, ...
																				fband, 'bandpass');
		%--------------------------------------------------
		% filter data using input filter settings
		%--------------------------------------------------
		resp = filtfilt(handles.fcoeffb, handles.fcoeffa, resp);
		%--------------------------------------------------
		% update analysis window
		%--------------------------------------------------
		% check if analysis window is beyond length of signal
		if  ms2samples(handles.Awindow(2), Fs) > length(resp)
			% if so, reset to duration of signal
			handles.Awindow(2) = floor(bin2ms(length(resp), Fs));
			update_ui_str(handles.AnalysisEndCtrl, handles.Awindow(2));
			guidata(hObject, handles);
			fprintf('warning: Analysis End > length of signal!!!!');
		end
		% find bins for analysis
		bin = ms2samples(handles.Awindow, Fs);
		if bin(1) == 0
			bin(1) = 1;
		end
		%--------------------------------------------------
		% analyze data
		%--------------------------------------------------
		% compute RMS
		resp_RMS = rms(resp(bin(1):bin(2)));
		% compute dB SPL
		resp_dBSPL = dbspl(handles.VtoPa*resp_RMS);
% 		% take fft of raw and adj response data
% 		[fresp, magresp, phiresp] = daqdbfullfft(	resp(bin(1):bin(2)), ...
% 																Fs, ...
% 																length(resp(bin(1):bin(2))) );
		%--------------------------------------------------
		% assign data to handles containers, update display
		%--------------------------------------------------
		if strcmpi(ButtonID, 'Play Raw')
			handles.rawresp = resp;
			handles.rawRMS = resp_RMS;
			handles.rawdBSPL = resp_dBSPL;
			if strcmpi(handles.PlotSpectrum, 'on')
				% take fft of raw response data
				[handles.rawfresp, handles.rawmag, handles.rawphi] = ...
											daqdbfullfft(	resp(bin(1):bin(2)), ...
																Fs, ...
																length(resp(bin(1):bin(2))) );
				% set limits for Mag and Phase plots
				dblim = [min(handles.rawmag) max(handles.rawmag)];
				freqlim = 0.001*[0 Fs/2];

			else
				handles.rawfresp = [];
				handles.rawmag = [];
				handles.rawphi = [];
			end
			% update display
			dbtext = sprintf('Raw dB SPL: %.2f  [%d - %d]\n', ...
																handles.rawdBSPL, ...
																handles.Awindow(1), ...
																handles.Awindow(2));
			fprintf('%s\n', dbtext);
			update_ui_str(handles.RawdBText, dbtext);
			show_uictrl(handles.RawdBText);
			guidata(hObject, handles);
		elseif strcmpi(ButtonID, 'Play Adj')
			handles.adjresp = resp;
			handles.adjRMS = resp_RMS;
			handles.adjdBSPL = resp_dBSPL;
			if strcmp(handles.PlotSpectrum, 'on')
				% take fft of response data
				[handles.adjfresp, handles.adjmag, handles.adjphi] = ...
											daqdbfullfft(	resp(bin(1):bin(2)), ...
																Fs, ...
																length(resp(bin(1):bin(2))) );
				% set limits for Mag and Phase plots
				dblim = [min(handles.adjmag) max(handles.adjmag)];
				freqlim = 0.001*[0 Fs/2];

			else
				handles.adjfresp = [];
				handles.adjmag = [];
				handles.adjphi = [];
			end
			% update display
			dbtext = sprintf('Adj dB SPL: %.2f  [%d - %d]\n', ...
																handles.adjdBSPL, ...
																handles.Awindow(1), ...
																handles.Awindow(2));
			fprintf('%s\n', dbtext);
			update_ui_str(handles.AdjdBText, dbtext);
			show_uictrl(handles.AdjdBText);
			guidata(hObject, handles);
		end
		guidata(hObject, handles);
		handles.Awindow
		%-----------------------------------------------------------------------
		% plot data
		%-----------------------------------------------------------------------
		if strcmpi(ButtonID, 'Play Raw')
			% raw plots
			subplot(handles.P.rsig)
			tvec = 1000 * (0:(length(resp)-1)) ./ Fs;
			plot(handles.P.rsig, tvec, resp)
			title(handles.P.rsig, 'Response (V)')
			ylabel(handles.P.rsig, 'Raw', 'Color', 'b')
			set(handles.P.rsig, 'XTickLabel', []);

			if strcmpi(handles.PlotSpectrum, 'on')
				subplot(handles.P.rmag)
				plot(handles.P.rmag, 0.001*handles.rawfresp, handles.magraw);
				title(handles.P.rmag, 'Magnitude (dB)')
				ylim(handles.P.rmag, dblim);
				xlim(handles.P.rmag, freqlim);
				set(handles.P.rmag, 'XTickLabel', []);

				subplot(handles.P.rphi)
				plot(handles.P.rphi, 0.001*handles.rawfresp, unwrap(handles.rawphi));
				title(handles.P.rphi, 'Phase (rad)')
				xlim(handles.P.rphi, freqlim);
				set(handles.P.rphi, 'XTickLabel', []);

				subplot(handles.P.rspec)
				[S, F, T, P] = spectrogram(	resp, ...
														handles.SpectrumWindow, ...
														floor(0.98*handles.SpectrumWindow), ...
														512, ...
														Fs	);
				surf(handles.P.rspec, 1000*T, 0.001*F, 20*log10(P), ...
							'edgecolor', 'none');
				ylim(handles.P.rspec, freqlim);
				axis(handles.P.rspec, 'tight');
				view(handles.P.rspec, 0, 90);
				title(handles.P.rspec, 'Time vs. Freq (kHz) vs. dB')
				set(handles.P.rspec, 'XTickLabel', []);
				colormap(handles.P.rspec, handles.ColorMap);
			else
				cla(handles.P.rmag);
				cla(handles.P.rphi);
				cla(handles.P.rspec);
			end
			guidata(hObject, handles);
			drawnow
			
			updateDBplots(hObject, eventdata, handles);
			if handles.dBPlot
				updateDBplots(hObject, eventdata, handles);
			end
			guidata(hObject, handles);

		elseif strcmpi(ButtonID, 'Play Adj')
			% Update adj plots
			subplot(handles.P.asig)
			tvec = 1000 * (0:(length(resp)-1)) ./ Fs;
			plot(handles.P.asig, tvec, resp, 'r')
			ylabel(handles.P.asig, 'Adj', 'Color', 'r')
			xlabel(handles.P.asig, 'time (ms)')

			if strcmpi(handles.PlotSpectrum, 'on')
				subplot(handles.P.amag)
				plot(handles.P.amag, 0.001*handles.adjfresp, handles.adjmag, 'r');
				ylim(handles.P.amag, dblim);
				xlim(handles.P.amag, freqlim);
				xlabel(handles.P.amag, 'freq (kHz)');

				subplot(handles.P.aphi)
				plot(handles.P.aphi, 0.001*handles.adjfresp, unwrap(handles.adjphi), 'r');
				xlim(handles.P.aphi, freqlim);
				xlabel(handles.P.aphi, 'freq (kHz)');

				subplot(handles.P.aspec)
				[S, F, T, P] = spectrogram(	resp, ...
														handles.SpectrumWindow, ...
														floor(0.98*handles.SpectrumWindow), ...
														512, ...
														Fs	);
				surf(handles.P.aspec, 1000*T, 0.001*F, 20*log10(P), 'edgecolor', 'none');
				ylim(handles.P.aspec, freqlim);
				axis(handles.P.aspec, 'tight')
				view(handles.P.aspec, 0, 90);
				xlabel(handles.P.aspec, 'Time (ms)')
				colormap(handles.P.aspec, handles.ColorMap);
			else
				cla(handles.P.amag);
				cla(handles.P.aphi);
				cla(handles.P.aspec);
			end
			guidata(hObject, handles);
			drawnow
			
			if handles.dBPlot
				dBPlotCtrl_Callback(hObject, eventdata, handles);
			end
			guidata(hObject, handles);
		end
		%------------------------------------------------------------------------
	end
	% assign outputs
	varargout{1} = hObject;
	varargout{2} = handles;
%------------------------------------------------------------------------------
%------------------------------------------------------------------------------


%------------------------------------------------------------------------------
%------------------------------------------------------------------------------
function dBPlotCtrl_Callback(hObject, eventdata, handles)
%------------------------------------------------------------------------------
% displays dB level along with response signal
	handles.dBPlot = read_ui_val(handles.dBPlotCtrl);
	guidata(hObject, handles);
	if handles.dBPlot
		[dBfigure, dBaxes] = updateDBplots(hObject, eventdata, handles);
		% store returned handles
		if ~isempty(dBfigure)
			handles.dBFigure = dBfigure;
		end
		if ~isempty(dBaxes.rawdb)
			handles.P.rawdb = dBaxes.rawdb;
		end
		if ~isempty(dBaxes.adjdb)
			handles.P.adjdb = dBaxes.adjdb;
		end
	end
	guidata(hObject, handles);
%------------------------------------------------------------------------------
%------------------------------------------------------------------------------

%------------------------------------------------------------------------------
%------------------------------------------------------------------------------
function [dBFigure, dBAxes] = updateDBplots(hObject, eventdata, handles)
%------------------------------------------------------------------------------
% updates dB plot
%
%	should probably decimate data plots!
%	
%------------------------------------------------------------------------
	%------------------------------------------------------------------------
	% check to make sure there are data in the handles...
	%------------------------------------------------------------------------
	if any(	[	isempty(handles.respFs) ...
					(isempty(handles.rawresp) && isempty(handles.adjresp))	])
		warning('FlatWav:Data', 'No data for db analysis');
		dBFigure = []; %#ok<NASGU>
		dBAxes = [];
	end
	%------------------------------------------------------------------------	
	% if dB figure has not been created, do so
	%------------------------------------------------------------------------
	if isempty(handles.dBFigure) || ~ishandle(handles.dBFigure)
		fprintf('Creating dBFigure\n\n');
		dBFigure = figure;
	else
		dBFigure = handles.dBFigure;
		figure(dBFigure);
	end
	%------------------------------------------------------------------------
	% if subplots don't exist, create them
	%------------------------------------------------------------------------
	if (~isfield(handles.P, 'rawdb') || ~isfield(handles.P, 'adjdb'))
		dBAxes.rawdb =	subplot(211);
		dBAxes.adjdb = subplot(212);
	elseif (~ishandle(handles.P.rawdb) || ~ishandle(handles.P.adjdb))
		dBAxes.rawdb =	subplot(211);
		dBAxes.adjdb = subplot(212);
	else
		dBAxes.rawdb = handles.P.rawdb;
		dBAxes.adjdb = handles.P.adjdb;
	end
	%------------------------------------------------------------------------
	% compute rmsbins window
	%------------------------------------------------------------------------
	if ~isempty(handles.respFs)
		rmsbins = ms2samples(handles.PeakRMSWindow, handles.respFs); %#ok<NASGU>
	else
		error('updateDBplots: handles.respFs is empty!\n\n');
	end
	%------------------------------------------------------------------------
	% compute rms of raw response, plot
	%------------------------------------------------------------------------
	if ~isempty(handles.rawresp)
		 plotSignalAnddB(	handles.rawresp, ...
								handles.PeakRMSWindow, ...
								handles.respFs, ...
								'Axes', dBAxes.rawdb, ...
								'dBSPL', handles.VtoPa, ...
								'SignalName', 'Raw', ...
								'SignalColor', 'b', ...
								'dBMarkerColor', 'b');		
	end
	%------------------------------------------------------------------------	
	% compute rms of adj response, plot
	%------------------------------------------------------------------------
	if ~isempty(handles.adjresp)
		
		 plotSignalAnddB(	handles.adjresp, ...
								handles.PeakRMSWindow, ...
								handles.respFs, ...
								'Axes', dBAxes.adjdb, ...
								'dBSPL', handles.VtoPa, ...
								'SignalName', 'Adj', ...
								'SignalColor', 'r', ...
								'dBMarkerColor', 'r');
	end
%------------------------------------------------------------------------------
%------------------------------------------------------------------------------

%------------------------------------------------------------------------------
%------------------------------------------------------------------------------
function varargout = plotSignalAnddB(signal, rmswin, Fs, varargin)
%------------------------------------------------------------------------------
	% definitions
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
	
	% check input arguments
	nvararg = length(varargin);
	if nvararg
		aindex = 1;
		while aindex <= nvararg
			switch(upper(varargin{aindex}))

				% select axes
				case 'AXES'
					if ~ishandle(varargin{aindex + 1})
						error('%s: invalid axes handle', mfilename);
					end
					dBAxes = varargin{aindex + 1};
					aindex = aindex + 2;

				% set dBSPL option
				case 'DBSPL'
					% store conversion factor
					VtoPa = varargin{aindex + 1};
					aindex = aindex + 2;

				% set signal trace name
				case 'SIGNALNAME'
					sigName = varargin{aindex + 1};
					aindex = aindex + 2;
					
				% set signal trace color
				case 'SIGNALCOLOR'
					sigColor = varargin{aindex + 1};
					aindex = aindex + 2;				
					
				% set signal trace style
				case 'SIGNALSTYLE'
					sigStyle = varargin{aindex + 1};
					aindex = aindex + 2;				

				% set dB trace color
				case 'DBCOLOR'
					dbColor = varargin{aindex + 1}; %#ok<NASGU>
					aindex = aindex + 2;				
					
				% set db trace style
				case 'DBSTYLE'
					dbStyle = varargin{aindex + 1}; %#ok<NASGU>
					aindex = aindex + 2;				
				
				% set db trace line width
				case 'DBLINEWIDTH'
					dBLineWidth = varargin{aindex + 1};
					aindex = aindex + 2;				
				
				% set db peak marker symbol
				case 'DBMARKER'
					dBMarker = varargin{aindex + 1};
					aindex = aindex + 2;				

				% set db peak marker size (points)
				case 'DBMARKERSIZE'
					dBMarkerSize = varargin{aindex + 1};
					aindex = aindex + 2;
				
				% set db peak marker size (points)
				case 'DBMARKERCOLOR'
					dBMarkerColor = varargin{aindex + 1};
					aindex = aindex + 2;				

				otherwise
					error('%s: Unknown option %s', mfilename, varargin{aindex});
			end		% END SWITCH
		end		% END WHILE aindex

	else
		dBFigure = figure;
		dBAxes = axes;	
	end		% END IF nvararg
	
	% comvert rmswindow from milliseconds to # of bins
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

	if nargout
		varargout{1} = dBFigure;
		varargout{2} = dBAxes;
	end
%------------------------------------------------------------------------------
%------------------------------------------------------------------------------
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
	[figfile, figpath] = ...
							uiputfile('*.fig','Save plot and figure in .fig file...');
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
		ylim(handles.CalibrationAxes, ...
				[0.9*min(handles.cal.mag(1, :)) 1.1*max(handles.cal.mag(1, :))]);
		grid(	handles.CalibrationAxes, 'on');
	end
	guidata(hObject, handles);
	SmoothCalCtrl_Callback(hObject, eventdata, handles);
%-------------------------------------------------------------------------

%-------------------------------------------------------------------------
function FlatCalMenuItem_Callback(hObject, eventdata, handles)
%--------------------------------------------------
% fake cal data
%--------------------------------------------------
	handles.cal = fake_caldata('freqs', (1:10:(handles.S.Fs / 2)));
	handles.cal.mag = 90 * handles.cal.mag;
	guidata(hObject, handles);
	plot(handles.CalibrationAxes, 0.001*handles.cal.freq, ...
																handles.cal.mag(1, :), '.-');
	ylim(handles.CalibrationAxes, [0 100]);
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
			save_audio(0.9*normalize(handles.adj), handles.S.Fs, datafile);
% 			peakfile = [datafile(1:(end-4)) '_PeakVal.txt'];
% 			save(peakfile, peakval, '-ascii');
		else
			save_audio(handles.adj, handles.S.Fs, datafile);
		end
	end
	guidata(hObject, handles);
%-------------------------------------------------------------------------

%-------------------------------------------------------------------------
function save_audio(audiodata, Fs, audiofile)
%-------------------------------------------------------------------------
	if verLessThan('matlab', 'R2015b')
		wavwrite(audiodata, Fs, audiofile); %#ok<DWVWR>
	else
		audiowrite(audiofile, audiodata, Fs);
	end
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
			save_audio(0.9*normalize(handles.raw), handles.S.Fs, datafile);
% 			peakfile = [datafile(1:(end-4)) '_PeakVal.txt'];
% 			save(peakfile, peakval, '-ascii');
		else
			save_audio(handles.raw, handles.S.Fs, datafile);
		end
	end
%-------------------------------------------------------------------------

%-------------------------------------------------------------------------
function SaveAllSignalsMenuItem_Callback(hObject, eventdata, handles)
%-------------------------------------------------------------------------
	[matfile, matpath] = uiputfile(	'*.mat', ...
												'Save signals to mat file...');
	if matfile ~= 0
		raw = handles.raw; %#ok<NASGU>
		adj = handles.adj; %#ok<NASGU>
		S = handles.S; %#ok<NASGU>
		if strcmpi(handles.SignalMode, 'WAV')
			wavdata = handles.wavdata; %#ok<NASGU>
			save(fullfile(matpath, matfile), 'raw', 'adj', 'S', 'wavdata', '-MAT');
			clear raw adj S wavdata;
		else
			synth = handles.synth; %#ok<NASGU>
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
	handles.HPFc = newVal;
	%{
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
	updatePlots(hObject, eventdata, handles);
%-------------------------------------------------------------------------


%-------------------------------------------------------------------------
function RMSWindowMenuItem_Callback(hObject, eventdata, handles)
	% get new spectrum window
	newVal = uiaskvalue(	...
					'Value',				handles.PeakRMSWindow,			...
					'ValueText',		'RMS window (ms)',		...
					'QuestionText',	'Input Window Size for Level Calculation', ...
					'FigureName',		''	);
	if ~isnumeric(newVal)
		warndlg('RMS size must be a number!', 'FlatWav')
	elseif ~between(newVal, 1, 1e6)
		warndlg('RMS Window size must be between 1 and 1e6!', 'FlatWav')
	else
		handles.PeakRMSWindow = newVal;
	end
	guidata(hObject, handles);
	updateDBplots(hObject, eventdata, handles)
	guidata(hObject, handles);
%-------------------------------------------------------------------------


%-------------------------------------------------------------------------
function PlotSpectrumMenuItem_Callback(hObject, eventdata, handles)
	oldval = get(handles.PlotSpectrumMenuItem, 'Checked');
	if strcmp(oldval, 'on')
		handles.PlotSpectrum = 'off';
	else
		handles.PlotSpectrum = 'on';
	end
	set(handles.PlotSpectrumMenuItem, 'Checked', handles.PlotSpectrum);
	guidata(hObject, handles);
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
function CompMethodCtrl_CreateFcn(hObject, eventdata, handles) %#ok<*INUSD>
	if ispc && isequal(get(hObject,'BackgroundColor'), ...
			get(0,'defaultUicontrolBackgroundColor'))
		  set(hObject,'BackgroundColor','white');
	end
	set(hObject, 'String', {'plot(rand(5))', 'plot(sin(1:0.01:25))', ...
		'bar(1:.5:10)', 'plot(membrane)', 'surf(peaks)'});
function WavFilenameCtrl_CreateFcn(hObject, eventdata, handles)
	if ispc && isequal(get(hObject,'BackgroundColor'), ...
			get(0,'defaultUicontrolBackgroundColor'))
		 set(hObject,'BackgroundColor','white');
	end
function SynthTypeCtrl_CreateFcn(hObject, eventdata, handles)
	if ispc && isequal(get(hObject,'BackgroundColor'), ...
			get(0,'defaultUicontrolBackgroundColor'))
		 set(hObject,'BackgroundColor','white');
	end
function FsCtrl_CreateFcn(hObject, eventdata, handles)
	if ispc && isequal(get(hObject,'BackgroundColor'), ...
			get(0,'defaultUicontrolBackgroundColor'))
		 set(hObject,'BackgroundColor','white');
	end
function p1Ctrl_CreateFcn(hObject, eventdata, handles)
	if ispc && isequal(get(hObject,'BackgroundColor'), ...
			get(0,'defaultUicontrolBackgroundColor'))
		 set(hObject,'BackgroundColor','white');
	end
function p2Ctrl_CreateFcn(hObject, eventdata, handles)
	if ispc && isequal(get(hObject,'BackgroundColor'), ...
			get(0,'defaultUicontrolBackgroundColor'))
		 set(hObject,'BackgroundColor','white');
	end
function p3Ctrl_CreateFcn(hObject, eventdata, handles)
	if ispc && isequal(get(hObject,'BackgroundColor'), ...
			get(0,'defaultUicontrolBackgroundColor'))
		 set(hObject,'BackgroundColor','white');
	end
function p4Ctrl_CreateFcn(hObject, eventdata, handles)
	if ispc && isequal(get(hObject,'BackgroundColor'), ...
			get(0,'defaultUicontrolBackgroundColor'))
		 set(hObject,'BackgroundColor','white');
	end
function p5Ctrl_CreateFcn(hObject, eventdata, handles)
	if ispc && isequal(get(hObject,'BackgroundColor'), ...
			get(0,'defaultUicontrolBackgroundColor'))
		 set(hObject,'BackgroundColor','white');
	end
function p6Ctrl_CreateFcn(hObject, eventdata, handles)
	if ispc && isequal(get(hObject,'BackgroundColor'), ...
			get(0,'defaultUicontrolBackgroundColor'))
		 set(hObject,'BackgroundColor','white');
	end
function TargetSPLCtrl_CreateFcn(hObject, eventdata, handles)
	if ispc && isequal(get(hObject,'BackgroundColor'), ...
			get(0,'defaultUicontrolBackgroundColor'))
		 set(hObject,'BackgroundColor','white');
	end
function LowCutFreqCtrl_CreateFcn(hObject, eventdata, handles)
	if ispc && isequal(get(hObject,'BackgroundColor'), ...
			get(0,'defaultUicontrolBackgroundColor'))
		 set(hObject,'BackgroundColor','white');
	end
function SpectrumWindowCtrl_CreateFcn(hObject, eventdata, handles)
	if ispc && isequal(get(hObject,'BackgroundColor'), ...
			get(0,'defaultUicontrolBackgroundColor'))
		 set(hObject,'BackgroundColor','white');
	end
function CorrFminCtrl_CreateFcn(hObject, eventdata, handles)
	if ispc && isequal(get(hObject,'BackgroundColor'), ...
			get(0,'defaultUicontrolBackgroundColor'))
		 set(hObject,'BackgroundColor','white');
	end
function CorrFmaxCtrl_CreateFcn(hObject, eventdata, handles)
	if ispc && isequal(get(hObject,'BackgroundColor'), ...
			get(0,'defaultUicontrolBackgroundColor'))
		 set(hObject,'BackgroundColor','white');
	end
function NormalizePeakCtrl_CreateFcn(hObject, eventdata, handles)
	if ispc && isequal(get(hObject,'BackgroundColor'), ...
			get(0,'defaultUicontrolBackgroundColor'))
		 set(hObject,'BackgroundColor','white');
	end
function AnalysisStartCtrl_CreateFcn(hObject, eventdata, handles)
	if ispc && isequal(get(hObject,'BackgroundColor'), ...
			get(0,'defaultUicontrolBackgroundColor'))
		 set(hObject,'BackgroundColor','white');
	end
function AnalysisEndCtrl_CreateFcn(hObject, eventdata, handles)
	if ispc && isequal(get(hObject,'BackgroundColor'), ...
			get(0,'defaultUicontrolBackgroundColor'))
		 set(hObject,'BackgroundColor','white');
	end
function SmoothVal1Ctrl_CreateFcn(hObject, eventdata, handles)
	if ispc && isequal(get(hObject,'BackgroundColor'), ...
			get(0,'defaultUicontrolBackgroundColor'))
		 set(hObject,'BackgroundColor','white');
	end
function CalSmoothMethodCtrl_CreateFcn(hObject, eventdata, handles)
	if ispc && isequal(get(hObject,'BackgroundColor'), ...
			get(0,'defaultUicontrolBackgroundColor'))
		set(hObject,'BackgroundColor','white');
	end
function SmoothVal2Ctrl_CreateFcn(hObject, eventdata, handles)
	if ispc && isequal(get(hObject,'BackgroundColor'), ...
			get(0,'defaultUicontrolBackgroundColor'))
	    set(hObject,'BackgroundColor','white');
	end
function PreFilterRangeCtrl_CreateFcn(hObject, eventdata, handles)
	if ispc && isequal(get(hObject,'BackgroundColor'), ...
			get(0,'defaultUicontrolBackgroundColor'))
		 set(hObject,'BackgroundColor','white');
	end
function PostFilterRangeCtrl_CreateFcn(hObject, eventdata, handles)
	if ispc && isequal(get(hObject,'BackgroundColor'), ...
			get(0,'defaultUicontrolBackgroundColor'))
		 set(hObject,'BackgroundColor','white');
	end
function CorrectionLimitValCtrl_CreateFcn(hObject, eventdata, handles)
	if ispc && isequal(get(hObject,'BackgroundColor'), ...
			get(0,'defaultUicontrolBackgroundColor'))
		 set(hObject,'BackgroundColor','white');
	end
function SmoothEdgesValCtrl_CreateFcn(hObject, eventdata, handles)
	if ispc && isequal(get(hObject,'BackgroundColor'), ...
			get(0,'defaultUicontrolBackgroundColor'))
		 set(hObject,'BackgroundColor','white');
	end
function PeakTimeRawCtrl_CreateFcn(hObject, eventdata, handles)
	if ispc && isequal(get(hObject,'BackgroundColor'), ...
			get(0,'defaultUicontrolBackgroundColor'))
		 set(hObject,'BackgroundColor','white');
	end
function PeakTimeAdjCtrl_CreateFcn(hObject, eventdata, handles)
	if ispc && isequal(get(hObject,'BackgroundColor'), ...
			get(0,'defaultUicontrolBackgroundColor'))
		 set(hObject,'BackgroundColor','white');
	end
function CalDisplayCtrl_CreateFcn(hObject, eventdata, handles)
	if ispc && isequal(get(hObject,'BackgroundColor'), ...
			get(0,'defaultUicontrolBackgroundColor'))
		 set(hObject,'BackgroundColor','white');
	end
	
%******************************************************************************
%******************************************************************************
%******************************************************************************

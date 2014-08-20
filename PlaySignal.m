function varargout = PlaySignal(hObject, handles)
%------------------------------------------------------------------------------
% PlaySignal
%------------------------------------------------------------------------------
% sets up NI data acquisition toolbox parameters and plays signal
%------------------------------------------------------------------------------

%------------------------------------------------------------------------------
% Sharad J Shanbhag
% sshanbhag@neomed.edu
%------------------------------------------------------------------------------
% Created: 4 December, 2012 (SJS)
% 				Created from NICal_NIinit.m and other scraps
% 
% Revisions:
%------------------------------------------------------------------------------

%--------------------------------------------------
% find who called us
%--------------------------------------------------
ButtonID = read_ui_str(hObject);

%--------------------------------------------------
% if figure has not been created, do so
%--------------------------------------------------


%--------------------------------------------------
% update bandpass filter for processing the data
%--------------------------------------------------
% Nyquist frequency
fnyq = handles.S.Fs / 2;
% passband definition
fband = [handles.HPFc handles.LPFc] ./ fnyq;
% filter coefficients using a butterworth bandpass filter
[handles.fcoeffb, handles.fcoeffa] = butter(handles.FilterOrder, ...
																		fband, 'bandpass');
% update analysis window
% check if analysis window is beyond length of signal
if  ms2samples(handles.Awindow(2), handles.S.Fs) > length(handles.raw)
	% if so, reset to duration of signal
	handles.Awindow(2) = floor(bin2ms(length(handles.raw), handles.S.Fs));
	% and update display
	update_ui_str(handles.AnalysisEndCtrl, handles.Awindow(2));
	guidata(hObject, handles);
	fprintf('warning: Analysis End > length of signal!!!!');
end
% find bins for analysis
bin = ms2samples(handles.Awindow, handles.S.Fs);
if bin(1) == 0
	bin(1) = 1;
end
% store changes to handles 
guidata(hObject, handles);

%--------------------------------------------------
% play sound out of selected output device
%--------------------------------------------------
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
	
elseif strcmpi(handles.OutputDevice, 'NIDAQ')
	if strcmpi(ButtonID, 'Play Raw') && ~isempty(handles.raw)
		% play raw sound
		disable_ui(hObject);
		handles.rawresp = NIplaysignal(hObject, handles);
		% compute dBSPL
		rawRMS = rms(handles.rawresp(bin(1):bin(2)));
		rawdBSPL = dbspl(handles.VtoPa*rawRMS);
		fprintf('Raw dB SPL: %.4f\n', rawdBSPL);
		update_ui_str(handles.RawdBText, sprintf('Raw dB SPL: %.2f', rawdBSPL));
		show_uictrl(handles.RawdBText);
		enable_ui(hObject);		
	elseif strcmpi(ButtonID, 'Play Adj') && ~isempty(handles.adj)
		% play adj sound
		disable_ui(hObject);
		handles.adjresp = NIplaysignal(hObject, handles);
		% compute dBSPL
		adjRMS = rms(handles.adjresp(bin(1):bin(2)));
		adjdBSPL = dbspl(handles.VtoPa*adjRMS);
		fprintf('Adj dB SPL: %.4f\n', adjdBSPL);
		update_ui_str(handles.AdjdBText, sprintf('Adj dB SPL: %.2f', adjdBSPL));
		show_uictrl(handles.AdjdBText);		
		enable_ui(hObject);
	end
	
else
	errordlg(sprintf('unknown io device %s', handles.OutputDevice), 'FlatWav Error');
end
guidata(hObject, handles);

% assign outputs
varargout{1} = hObject;
varargout{2} = handles;

return


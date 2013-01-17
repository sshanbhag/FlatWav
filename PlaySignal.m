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
ButtonID = read_ui_str(hObject) %#ok<NOPRT>

%--------------------------------------------------
% update bandpass filter for processing the data
%--------------------------------------------------
% Nyquist frequency
fnyq = handles.S.Fs / 2;
% passband definition
fband = [handles.HPFc handles.LPFc] ./ fnyq;
% filter coefficients using a butterworth bandpass filter
[handles.fcoeffb, handles.fcoeffa] = butter(handles.FilterOrder, fband, 'bandpass');
guidata(hObject, handles);

if strcmpi(handles.OutputDevice, 'WINSOUND')
	if (handles.S.Fs == 44100)  
		if strcmpi(ButtonID, 'Play Raw') && ~isempty(handles.raw)
			disable_ui(hObject);
			update_ui_str(hObject, 'RAW');
			pause(0.5);
			sound(sin2array(handles.raw, 1, handles.S.Fs), handles.S.Fs);
			pause(1);
			enable_ui(hObject);
			update_ui_str(handles.RawdBText, sprintf('Raw dB SPL: ---'));
			hide_uictrl(handles.RawdBText);
			hide_uictrl(handles.AdjdBText);
		elseif strcmpi(ButtonID, 'Play Adj') && ~isempty(handles.adj)
			disable_ui(hObject);
			pause(0.5);
			sound(sin2array(handles.adj, 1, handles.S.Fs), handles.S.Fs);
			pause(1);
			enable_ui(hObject);
			update_ui_str(handles.AdjdBText, sprintf('Adj dB SPL: ---'));
			hide_uictrl(handles.RawdBText);
			hide_uictrl(handles.AdjdBText);		
		end
	end
elseif strcmpi(handles.OutputDevice, 'NIDAQ')
	if strcmpi(ButtonID, 'Play Raw') && ~isempty(handles.raw)	
		handles.rawresp = NIplaysignal(hObject, handles);
		% compute dBSPL
		rawRMS = rms(handles.rawresp);
		rawdBSPL = dbspl(handles.VtoPa*rawRMS);
		fprintf('Raw dB SPL: %.4f\n', rawdBSPL);
		update_ui_str(handles.RawdBText, sprintf('Raw dB SPL: %.2f', rawdBSPL));
		show_uictrl(handles.RawdBText);
	elseif strcmpi(ButtonID, 'Play Adj') && ~isempty(handles.adj)
		handles.adjresp = NIplaysignal(hObject, handles);
		% compute dBSPL
		adjRMS = rms(handles.adjresp);
		adjdBSPL = dbspl(handles.VtoPa*adjRMS);
		fprintf('Adj dB SPL: %.4f\n', adjdBSPL);
		update_ui_str(handles.AdjdBText, sprintf('Adj dB SPL: %.2f', adjdBSPL));
		show_uictrl(handles.AdjdBText);		
	end

else
	errordlg(sprintf('unknown io device %s', handles.OutputDevice), 'FlatWav Error');
end
guidata(hObject, handles);

varargout{1} = hObject;
varargout{2} = handles;


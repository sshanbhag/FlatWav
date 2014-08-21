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
% setup output figure
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
	pw = 0.19;
	ph = 0.375;
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
	% take fft of raw and adj response data
	[fresp, magresp, phiresp] = daqdbfullfft(	resp(bin(1):bin(2)), ...
															Fs, ...
															length(resp(bin(1):bin(2))) );
	
	%--------------------------------------------------
	% assign data to handles containers, update display
	%--------------------------------------------------
	if strcmpi(ButtonID, 'Play Raw')
		handles.rawresp = resp;
		handles.rawRMS = resp_RMS;
		handles.rawdBSPL = resp_dBSPL;
		handles.rawfresp = fresp;
		handles.rawmag = magresp;
		handles.rawphi = phiresp;
		fprintf('Raw dB SPL: %.4f\n', handles.rawdBSPL);
		update_ui_str(handles.RawdBText, ...
							sprintf('Raw dB SPL: %.2f', handles.rawdBSPL));
		show_uictrl(handles.RawdBText);

	elseif strcmpi(ButtonID, 'Play Adj')
		handles.adjresp = resp;
		handles.adjRMS = resp_RMS;
		handles.adjdBSPL = resp_dBSPL;
		handles.adjfresp = fresp;
		handles.adjmag = magresp;
		handles.adjphi = phiresp;
		fprintf('Adj dB SPL: %.4f\n', handles.adjdBSPL);
		update_ui_str(handles.AdjdBText, ...
							sprintf('Adj dB SPL: %.2f', handles.adjdBSPL));
		show_uictrl(handles.AdjdBText);		
	end
	guidata(hObject, handles);
	
	%-----------------------------------------------------------------------
	% plot data
	%-----------------------------------------------------------------------
	% set limits for Mag and Phase plots
	dblim = [min(magresp) max(magresp)];
	freqlim = 0.001*[0 Fs/2];

	if strcmpi(ButtonID, 'Play Raw')
		% raw plots
		subplot(handles.P.rsig)
		tvec = 1000 * (0:(length(resp)-1)) ./ Fs;
		plot(tvec, resp)
		title('Response (V)')
		ylabel('Raw', 'Color', 'b')
		set(gca, 'XTickLabel', []);

		subplot(handles.P.rmag)
		plot(0.001*fresp, magresp);
		title('Magnitude (dB)')
		ylim(dblim);
		xlim(freqlim);
		set(gca, 'XTickLabel', []);

		subplot(handles.P.rphi)
		plot(0.001*fresp, unwrap(phiresp));
		title('Phase (rad)')
		xlim(freqlim);
		set(gca,  'XTickLabel', []);

		subplot(handles.P.rspec)
		[S, F, T, P] = spectrogram(	resp, ...
												handles.SpectrumWindow, ...
												floor(0.95*handles.SpectrumWindow), ...
												512, ...
												Fs	);
		surf(1000*T, 0.001*F, 20*log10(P), 'edgecolor', 'none');
		ylim(freqlim);
		axis tight;
		view(0, 90);
		title('Time vs. Freq (kHz) vs. dB')
		set(gca,  'XTickLabel', []);
		colormap(handles.P.rspec, handles.ColorMap)

	elseif strcmpi(ButtonID, 'Play Adj')
		% Update adj plots
		subplot(handles.P.asig)
		tvec = 1000 * (0:(length(resp)-1)) ./ Fs;
		plot(tvec, resp, 'r')
		ylabel('Adj', 'Color', 'r')
		xlabel('time (ms)')

		subplot(handles.P.amag)
		plot(0.001*fresp, magresp, 'r');
		ylim(dblim);
		xlim(freqlim);
		xlabel('freq (kHz)');

		subplot(handles.P.aphi)
		plot(0.001*fresp, unwrap(phiresp), 'r');
		xlim(freqlim);
		xlabel('freq (kHz)');

		subplot(handles.P.aspec)
		[S, F, T, P] = spectrogram(	resp, ...
												handles.SpectrumWindow, ...
												floor(0.95*handles.SpectrumWindow), ...
												512, ...
												Fs	);
		surf(1000*T, 0.001*F, 20*log10(P), 'edgecolor', 'none');
		ylim(freqlim);
		axis tight;
		view(0, 90);
		xlabel('Time (ms)')
		colormap(handles.P.aspec, handles.ColorMap);
	end
	%------------------------------------------------------------------------------
end

% assign outputs
varargout{1} = hObject;
varargout{2} = handles;

return


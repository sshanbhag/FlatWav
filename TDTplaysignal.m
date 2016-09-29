function varargout = TDTplaysignal(hObject, handles, outputsignal)
%------------------------------------------------------------------------------
% [resp, Fs] = TDTplaysignal(hObject, handles, outputsignal)
%------------------------------------------------------------------------------
% FlatWav
%------------------------------------------------------------------------------
% sets up TDT parameters and plays signal
%------------------------------------------------------------------------------

%------------------------------------------------------------------------------
% Sharad J Shanbhag
% sshanbhag@neomed.edu
%------------------------------------------------------------------------------
% Created: 29 September, 2016 (SJS)
% 				Created from NIplaysignal.m and other scraps
% 
% Revisions:
%------------------------------------------------------------------------------

%-----------------------------------------------------------------------
%-----------------------------------------------------------------------
% Settings/Constants
%-----------------------------------------------------------------------
%-----------------------------------------------------------------------
% local copy of TDT struct and iodev struct
TDT = handles.TDT;
iodev = handles.TDT.iodev;

%--------------------------------------------------
% limits to analog input and output
%--------------------------------------------------
% TDT_Constants;
AI_LIMIT = 10;
AO_LIMIT = 10;

%-----------------------------------------------------------------------
%-----------------------------------------------------------------------
% Check Inputs
%-----------------------------------------------------------------------
%-----------------------------------------------------------------------
% check to make sure output signal isn't crazy
if isempty(outputsignal)	
	warndlg('outputsignal is empty!');
	return
elseif max(abs(outputsignal)) > AO_LIMIT
	warndlg('outputsignal out of range');
	return
end

%-----------------------------------------------------------------------
%-----------------------------------------------------------------------
% Initialize the TDT device
%-----------------------------------------------------------------------
%-----------------------------------------------------------------------
fprintf('%s: settings for TDT hardware...\n', mfilename);

inpts = length(outputsignal) + ms2samples(10, iodev.Fs);
set(iodev.NI.ai, 'SamplesPerTrigger', inpts);

% set the total sweep period time
RPsettag(iodev, 'SwPeriod', inpts);
% set the sweep count (may not be necessary)
RPsettag(iodev, 'SwCount', 1);
% Set the length of time to acquire data
RPsettag(iodev, 'AcqDur', inpts);
% set the stimulus delay
RPsettag(iodev, 'StimDelay', 0);
% set the stimulus Duration
RPsettag(iodev, 'StimDur', length(outputsignal));

%------------------------------------------------------------------------
%------------------------------------------------------------------------
%% I/O 
%------------------------------------------------------------------------
%------------------------------------------------------------------------
atten_val = [20 20];
% set attenuators
if strcmpi(TDT.AttenMode, 'PA5')
	% no need to test attenuation but do need to set the attenuators
	TDT.setattenFunc([PA5L PA5R], atten_val);
elseif strcmpi(TDT.AttenMode, 'RZ6')
	TDT.setattenFunc(iodev, atten_val);
elseif strcmpi(TDT.AttenMode, 'DIGITAL')
	% do something
end

%-------------------------------------------------------
% send data to hardware and play it while recording response
%-------------------------------------------------------
% play the sound;
[resp, rate] = TDT.iofunc(iodev, outputsignal, inpts);
% pause
pause(handles.IOpause);

%-----------------------------------------------------------------------
%-----------------------------------------------------------------------
% done, return data
%-----------------------------------------------------------------------
%-----------------------------------------------------------------------
if any(nargout == 1:2)
	varargout{1} = resp;
end
if nargout == 2
 	varargout{2} = rate;
end


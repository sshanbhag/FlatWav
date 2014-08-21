function varargout = NIplaysignal(hObject, handles, outputsignal)
%------------------------------------------------------------------------------
% [resp, magresp, phiresp] = NIplaysignal(hObject, handles, outputsignal)
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
% 21 Aug 2014 (SJS): 
%	-	pulled out plotting of signals, moved to PlaySignal.m
% 		This required addition of outputsig input arg and modification
%		of output args
%------------------------------------------------------------------------------

%-----------------------------------------------------------------------
%-----------------------------------------------------------------------
% Settings/Constants
%-----------------------------------------------------------------------
%-----------------------------------------------------------------------

%--------------------------------------------------
% limits to analog input and output
%--------------------------------------------------
% NICal_Constants;
AI_LIMIT = 5;
AO_LIMIT = 10;
% range needs to be in [RangeMin RangeMax] format
aiRange = AI_LIMIT * [-1 1];
aoRange = AO_LIMIT * [-1 1];

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
% Initialize the NI device
%-----------------------------------------------------------------------
%-----------------------------------------------------------------------
fprintf('%s: starting NI hardware...\n', mfilename);
iodev.Dnum = 'Dev1';
try
	iodev.NI = nidaq_aiao_singlechannel_init('NI', iodev.Dnum);
catch errMsg
	disp('error initializing NI device')
	disp(errMsg.identifier)
	disp(errMsg.message)
	disp(errMsg.cause)
	disp(errMsg.stack)
	return
end

%-----------------------------------------------------------------------
%-----------------------------------------------------------------------
%% set sample rate to value specified in cal settings
%-----------------------------------------------------------------------
%-----------------------------------------------------------------------
% AI subsystem
set(iodev.NI.ai, 'SampleRate', handles.S.Fs);
ActualRate = get(iodev.NI.ai, 'SampleRate');
if handles.S.Fs ~= ActualRate
	warning('FlatWav:NIDAQ', 'Requested ai Fs (%f) ~= ActualRate (%f)', handles.S.Fs, ActualRate);
end
iodev.Fs = ActualRate;
handles.S.Fs = iodev.Fs;
guidata(hObject, handles);

% AO subsystem
set(iodev.NI.ao, 'SampleRate', iodev.Fs);
ActualRate = get(iodev.NI.ao, 'SampleRate');
if iodev.Fs ~= ActualRate
	warning('NICAl:NIDAQ', 'ao: Requested SampleRate (%f) ~= ActualRate (%f)', iodev.Fs, ActualRate);
end
iodev.Fs = ActualRate;

%-----------------------------------------------------------------------
%-----------------------------------------------------------------------
% set input range
%-----------------------------------------------------------------------
%-----------------------------------------------------------------------
% set analog input range (might be overkill to set 
% InputRange, SensorRange and UnitsRange, but is seems to work)
for n = 1:length(iodev.NI.ai.Channel)
	iodev.NI.ai.Channel(n).InputRange = aiRange;
	iodev.NI.ai.Channel(n).SensorRange = aiRange;
	iodev.NI.ai.Channel(n).UnitsRange = aiRange;
end
% set analog output range
for n = 1:length(iodev.NI.ao.Channel)
	iodev.NI.ao.Channel(n).OutputRange = aoRange;
	iodev.NI.ao.Channel(n).UnitsRange = aoRange;
end

%------------------------------------------------------------------------
%------------------------------------------------------------------------
% HARDWARE TRIGGERING
%------------------------------------------------------------------------
%------------------------------------------------------------------------
% set TriggerType to manual (to synchronize ai and ao)
set([iodev.NI.ai iodev.NI.ao], 'TriggerType', 'Manual');
% set manual trigger to HW on
set(iodev.NI.ai,'ManualTriggerHwOn','Trigger')
% only 1 "sweep" per trigger event 
set(iodev.NI.ai, 'TriggerRepeat', 0);
% set SamplesPerTrigger to Inf for continous acquisition or 
% to # of samples to collect for each trigger event
inpts = length(outputsignal) + ms2samples(10, iodev.Fs);
set(iodev.NI.ai, 'SamplesPerTrigger', inpts);

%------------------------------------------------------------------------
%------------------------------------------------------------------------
%% set logging mode
%	'Disk'	sets logging mode to a file on disk (specified by 'LogFileName)
%	'Memory'	sets logging mode to memory only
%	'Disk&Memory'	logs to file and memory
%------------------------------------------------------------------------
%------------------------------------------------------------------------
set(iodev.NI.ai, 'LoggingMode', 'Memory');

%------------------------------------------------------------------------
%------------------------------------------------------------------------
%% set channel skew mode to Equisample
%------------------------------------------------------------------------
%------------------------------------------------------------------------
set(iodev.NI.ai, 'ChannelSkewMode', 'Equisample');

%------------------------------------------------------------------------
%------------------------------------------------------------------------
%% I/O 
%------------------------------------------------------------------------
%------------------------------------------------------------------------

%-------------------------------------------------------
% send ramped data to hardware and play it
%-------------------------------------------------------
putdata(iodev.NI.ao, sin2array(outputsignal, 1, handles.S.Fs)');
% wait for time to settle
timeToWait = ceil(bin2seconds(inpts, iodev.Fs)*2);

% start
start([iodev.NI.ai iodev.NI.ao]);
trigger([iodev.NI.ai iodev.NI.ao]);
wait(iodev.NI.ai, timeToWait);
% stop acquiring
stop([iodev.NI.ai iodev.NI.ao]);
% read data from ai object
index = iodev.NI.ai.SamplesAvailable;
resp = getdata(iodev.NI.ai, index);
% pause
pause(handles.IOpause);

%-----------------------------------------------------------------------
%-----------------------------------------------------------------------
% Clean up the NI Device
%-----------------------------------------------------------------------
%-----------------------------------------------------------------------
disp('...closing NI devices...');

% get event log
EventLogAI = showdaqevents(iodev.NI.ai);
EventLogAO = showdaqevents(iodev.NI.ao);

% delete and clear ai and ch0 object
delete(iodev.NI.ai);
delete(iodev.NI.ao);
delete(iodev.NI.chI);
delete(iodev.NI.chO);
clear iodev.NI.ai iodev.NI.ao iodev.NI.chI iodev.NI.chO

% save log information to mat file
[folder, mname, mext] = fileparts(which('FlatWav'));
save(fullfile(folder, 'EventLogs.mat'), ...
		'EventLogAI'			, ...
		'EventLogAO'			, ...
		'-MAT' );

clear iodev

if any(nargout == 1:2)
	varargout{1} = resp;
end
if nargout == 2
 	varargout{2} = ActualRate;
end


function varargout = NIplaysignal(handles)
%------------------------------------------------------------------------------
% NIplaysignal
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
%% Settings/Constants
%-----------------------------------------------------------------------
%-----------------------------------------------------------------------
% NICal_Constants;
AI_LIMIT = 5;
AO_LIMIT = 10;

% range needs to be in [RangeMin RangeMax] format
aiRange = AI_LIMIT * [-1 1];
aoRange = AO_LIMIT * [-1 1];

%-----------------------------------------------------------------------
%-----------------------------------------------------------------------
%% Check Inputs
%-----------------------------------------------------------------------
%-----------------------------------------------------------------------
% check to make sure output signal isn't crazy
if max(abs(handles.raw)) > AO_LIMIT
	warndlg('RAW signal out of range');
	return
end
if max(abs(handles.adj)) > AO_LIMIT
	warndlg('ADJ signal out of range!');
	return
end


%-----------------------------------------------------------------------
%-----------------------------------------------------------------------
%% Initialize the NI device
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
% AO subsystem
set(iodev.NI.ao, 'SampleRate', iodev.Fs);
ActualRate = get(iodev.NI.ao, 'SampleRate');
if iodev.Fs ~= ActualRate
	warning('NICAl:NIDAQ', 'ao: Requested SampleRate (%f) ~= ActualRate (%f)', iodev.Fs, ActualRate);
end
iodev.Fs = ActualRate;

%-----------------------------------------------------------------------
%-----------------------------------------------------------------------
%% set input range
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
%% HARDWARE TRIGGERING
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
inpts = length(handles.raw) + ms2samples(10, iodev.Fs);
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
% send raw data to hardware and play it
%-------------------------------------------------------
putdata(iodev.NI.ao, sin2array(handles.raw, 1, handles.S.Fs)');
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
rawresp = getdata(iodev.NI.ai, index);

% pause
pause(1);

%-------------------------------------------------------
% send adj data to hardware and play it
%-------------------------------------------------------
putdata(iodev.NI.ao, sin2array(handles.adj, 1, handles.S.Fs)');
% wait for time to settle
timeToWait = bin2seconds(inpts, iodev.Fs)*2;
% start
start([iodev.NI.ai iodev.NI.ao]);
trigger([iodev.NI.ai iodev.NI.ao]);
wait(iodev.NI.ai, 5);
% stop acquiring
stop([iodev.NI.ai iodev.NI.ao]);
% read data from ai object
index = iodev.NI.ai.SamplesAvailable;
adjresp = getdata(iodev.NI.ai, index);

%-------------------------------------------------------
%% filter data
%-------------------------------------------------------
rawresp = filtfilt(handles.fcoeffb, handles.fcoeffa, rawresp);
adjresp = filtfilt(handles.fcoeffb, handles.fcoeffa, adjresp);

%-----------------------------------------------------------------------
%% plot data
%-----------------------------------------------------------------------
% take fft of raw and adj response data
[fraw, magraw, phiraw] = daqdbfullfft(rawresp, iodev.Fs, length(rawresp));
[fadj, magadj, phiadj] = daqdbfullfft(adjresp, iodev.Fs, length(adjresp));


% plotting limits
dblim = [-120 0];
freqlim = 0.001*[0 iodev.Fs/2];

% raw plots
axes(handles.RawSignalAxes)
tvec = 1000 * (0:(length(rawresp)-1)) ./ iodev.Fs;
plot(tvec, rawresp)
title('Signal (V)')
ylabel('Raw', 'Color', 'b')
set(handles.RawSignalAxes, 'XTickLabel', []);

axes(handles.RawMagAxes)
plot(0.001*fraw, magraw);
title('Magnitude (dB)')
ylim(dblim);
xlim(freqlim);
set(handles.RawMagAxes, 'XTickLabel', []);

axes(handles.RawPhaseAxes)
plot(0.001*fraw, unwrap(phiraw));
title('Phase (rad)')
xlim(freqlim);
set(handles.RawPhaseAxes, 'XTickLabel', []);

axes(handles.RawSpectrumAxes)
[S, F, T, P] = spectrogram(	rawresp, ...
										handles.SpectrumWindow, ...
										floor(0.95*handles.SpectrumWindow), ...
										512, ...
										iodev.Fs	);
surf(1000*T, 0.001*F, 20*log10(P), 'edgecolor', 'none');
ylim(freqlim);
axis tight;
view(0, 90);
title('Time vs. Freq (kHz) vs. dB')
set(handles.RawSpectrumAxes, 'XTickLabel', []);
colormap(handles.AdjSpectrumAxes, handles.ColorMap)

% Update adj plots
axes(handles.AdjSignalAxes)
tvec = 1000 * (0:(length(adjresp)-1)) ./ iodev.Fs;
plot(tvec, adjresp, 'g')
ylabel('Adj', 'Color', 'g')
xlabel('time (ms)')

axes(handles.AdjMagAxes)
plot(0.001*fadj, magadj, 'g');
ylim(dblim);
xlim(freqlim);
xlabel('freq (kHz)');

axes(handles.AdjPhaseAxes)
plot(0.001*fadj, unwrap(phiadj), 'g');
xlim(freqlim);
xlabel('freq (kHz)');

axes(handles.AdjSpectrumAxes)
[S, F, T, P] = spectrogram(	adjresp, ...
										handles.SpectrumWindow, ...
										floor(0.95*handles.SpectrumWindow), ...
										512, ...
										iodev.Fs	);
surf(1000*T, 0.001*F, 20*log10(P), 'edgecolor', 'none');
ylim(freqlim);
axis tight;
view(0, 90);
xlabel('Time (ms)')
colormap(handles.AdjSpectrumAxes, handles.ColorMap)
%------------------------------------------------------------------------------


%-----------------------------------------------------------------------
%% Clean up the NI Device
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

% save settings information to mat file
save(fullfile(pwd, 'EventLogs.mat'), ...
		'EventLogAI'			, ...
		'EventLogAO'			, ...
		'-MAT' );

clear iodev

nargout
if nargout <= 2
	varargout{1} = rawresp;
end
if nargout == 2
	varargout{2} = adjresp;
end


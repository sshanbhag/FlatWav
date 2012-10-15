function [iodev, init_status] = ni_ioinit(Dnum, SweepDuration, Fs, IOrange)
%--------------------------------------------------------------------------
% ni_ioinit.m
%--------------------------------------------------------------------------
% sets up NI data acquisition toolbox parameters
%
%	Dnum		device id (usually 'Dev1')
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Sharad J Shanbhag
% sshanbhag@neomed.edu
%--------------------------------------------------------------------------
% Created: 15 October 2012 (SJS)
% 				Created from NICal_NIinit.m
% 
% Revisions:
%--------------------------------------------------------------------------

fprintf('%s: starting NI hardware...\n', mfilename);

%-----------------------------------------------------------------------
%-----------------------------------------------------------------------
% Settings/Constants
%-----------------------------------------------------------------------
%-----------------------------------------------------------------------
init_status = 0;

% initialize iodev struct
iodev = struct(	'Dnum', Dnum, ...
						'NI', [], ...
						'init_status', init_status, ...
						'Fs', [] ...
					);

%-----------------------------------------------------------------------
%-----------------------------------------------------------------------
% Initialize the NI device
%-----------------------------------------------------------------------
%-----------------------------------------------------------------------
try
	iodev.NI = nidaq_aiao_init('NI', Dnum);
catch errMsg
	disp('error initializing NI device')
	disp(errMsg.identifier);
	return
end

%-----------------------------------------------------------------------
%-----------------------------------------------------------------------
% set sample rate to value specified in cal settings
%-----------------------------------------------------------------------
%-----------------------------------------------------------------------

%------------------------------------------------------
% AI subsystem
%------------------------------------------------------
set(iodev.NI.ai, 'SampleRate', Fs);
ActualRate = get(iodev.NI.ai, 'SampleRate');
if Fs ~= ActualRate
	warning('ni_ioinit:NIDAQ', 'Requested ai Fs (%f) ~= ActualRate (%f)', Fs, ActualRate);
end
iodev.Fs = ActualRate;

%------------------------------------------------------
% AO subsystem
%------------------------------------------------------
set(iodev.NI.ao, 'SampleRate', Fs);
ActualRate = get(iodev.NI.ao, 'SampleRate');
if Fs ~= ActualRate
	warning('ni_ioinit:NIDAQ', 'ao: Requested SampleRate (%f) ~= ActualRate (%f)', Fs, ActualRate);
elseif iodev.Fs ~= ActualRate
	warning('ni_ioinit:NIDAQ', 'ao: analog in rate (%f) ~= ActualRate (%f)', iodev.Fs, ActualRate);
end
iodev.Fs = ActualRate;

%-----------------------------------------------------------------------
%-----------------------------------------------------------------------
% set input range
%-----------------------------------------------------------------------
%-----------------------------------------------------------------------
% range needs to be in [RangeMin RangeMax] format
aiaoRange = IOrange * [-1 1];
% set analog input range (might be overkill to set 
% InputRange, SensorRange and UnitsRange, but is seems to work)
for n = 1:length(iodev.NI.ai.Channel)
	iodev.NI.ai.Channel(n).InputRange = aiaoRange;
	iodev.NI.ai.Channel(n).SensorRange = aiaoRange;
	iodev.NI.ai.Channel(n).UnitsRange = aiaoRange;
end
% set analog output range
for n = 1:length(iodev.NI.ao.Channel)
	iodev.NI.ao.Channel(n).OutputRange = aiaoRange;
	iodev.NI.ao.Channel(n).UnitsRange = aiaoRange;
end

%------------------------------------------------------------------------
% HARDWARE TRIGGERING
%------------------------------------------------------------------------
% set TriggerType to manual (to synchronize ai and ao)
set([iodev.NI.ai iodev.NI.ao], 'TriggerType', 'Manual');
% set manual trigger to HW on
set(iodev.NI.ai,'ManualTriggerHwOn','Trigger')
% only 1 "sweep" per trigger event 
set(iodev.NI.ai, 'TriggerRepeat', 0);
% set SamplesPerTrigger to Inf for continous acquisition or 
% to # of samples to collect for each trigger event
set(iodev.NI.ai, 'SamplesPerTrigger', ms2samples(SweepDuration, iodev.Fs));

%-------------------------------------------------------
% set logging mode
%	'Disk'	sets logging mode to a file on disk (specified by 'LogFileName)
%	'Memory'	sets logging mode to memory only
%	'Disk&Memory'	logs to file and memory
%-------------------------------------------------------
set(iodev.NI.ai, 'LoggingMode', 'Memory');

%-------------------------------------------------------
% set channel skew mode to Equisample
%-------------------------------------------------------
set(iodev.NI.ai, 'ChannelSkewMode', 'Equisample');

%-------------------------------------------------------
% set init_status to 1
%-------------------------------------------------------
init_status = 1;
iodev.init_status = init_status;




function FlatWav_TDTsettings_RZ6(iodev, cal)
%------------------------------------------------------------------------
% FlatWav_TDTsettings_RZ6(iodev, cal)
%------------------------------------------------------------------------
% sets up TDT tag settings for FlatWav using RZ6
% 
%------------------------------------------------------------------------
% Input Arguments:
% 	iodev			TDT device interface structure
% 	cal             calibration data structure
% 
% Output Arguments:
%
%------------------------------------------------------------------------

%------------------------------------------------------------------------
%  Sharad Shanbhag & Go Ashida
%	sshanbhag@neomed.edu
%   ashida@umd.edu
%------------------------------------------------------------------------
%------------------------------------------------------------------------
% Originally Written (HeadphoneCal_tdtinit): 2009-2011 by SJS
% Upgraded Version Written (HeadphoneCal2_TDTsettings): 2011-2012 by GA
% 
% Revisions:
% 28 September (SJS): HeadphoneCal2_TDT_settings_RZ6 copied to FlatWav
%------------------------------------------------------------------------

%npts = 150000;  % size of the serial buffer -- fixed
%mclock = config.RPgettagFunc(iodev, 'mClock');

% set the TTL pulse duration
% RPsettag(iodev, 'TTLPulseDur', ms2samples(cal.TTLPulseDur, iodev.Fs));
% set the total sweep period time
RPsettag(iodev, 'SwPeriod', ms2samples(cal.SweepPeriod, iodev.Fs));
% set the sweep count (may not be necessary)
RPsettag(iodev, 'SwCount', 1);
% Set the length of time to acquire data
RPsettag(iodev, 'AcqDur', ms2samples(cal.AcqDuration, iodev.Fs));
% set the stimulus delay
RPsettag(iodev, 'StimDelay', ms2samples(cal.Delay, iodev.Fs));
% set the stimulus Duration
RPsettag(iodev, 'StimDur', ms2samples(cal.Duration, iodev.Fs));


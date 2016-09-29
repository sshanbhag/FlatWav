function [outhandles, outflag] = TDTclose(config, iodev, PA5L, PA5R)
%------------------------------------------------------------------------
% [ outhandles, outflag ] = TDTclose(config, indev, outdev, zBUS, PA5L, PA5R)
%------------------------------------------------------------------------
% 
% Closes/shuts down TDT I/O Hardware for FlatWav program
% 
%------------------------------------------------------------------------
% Input Arguments:
%   config        matlab struct containing configuration information
%   iodev       matlab struct containing io device information
%   PA5L        matlab struct containing PA5 (left) device information
%   PA5R        matlab struct containing PA5 (right) device information
% Output Arguments:
%   outhandles  handle containing iodev, PA5L, PA5R
%   outflag    flag to show if TDT is successfully terminated 
%               -1: error
%                0: not terminated 
%                1: success
%------------------------------------------------------------------------

%------------------------------------------------------------------------
% Sharad Shanbhag
%   sshanbhag@neomed.edu
%------------------------------------------------------------------------
% Created 29 September, 2016 from opto_TDTclose.m
%------------------------------------------------------------------------

disp([ mfilename ': ...closing TDT devices...']);
outflag = 0; %#ok<NASGU> % not terminated

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% check the TDT lock file 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist(config.TDTLOCKFILE, 'file')
    disp([ mfilename ': TDT lock file not found: ', config.TDTLOCKFILE ]);
    disp('Creating lock file, assuming TDT hardware is not initialized');
    TDTINIT = 0;
    save(config.TDTLOCKFILE, 'TDTINIT');
else
    load(config.TDTLOCKFILE);  % load the lock information
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Exit gracefully (close TDT objects, etc)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if TDTINIT
	outhandles = struct(); %#ok<UNRCH>
	outhandles.iodev = indev;
	outhandles.PA5L = PA5L;
	outhandles.PA5R = PA5R;

    %------------------------------------------------------------------
    % setting zBUS/indev/outdev/PA5 structure accoring to the config info
    % setiodev() are defined below
    %------------------------------------------------------------------
    [ iodev] = setiodev;

    %------------------------------------------------------------------
    % terminate PA5, RX*, zBUS 
    % pa5.closeFunc, indev.closeFunc, zbus.closeFunc etc. are defined below
    %------------------------------------------------------------------
	 if ~isempty(outhandles.PA5L)
	    disp('...closing PA5L')
	    outhandles.PA5L.status = pa5.closeFunc(PA5L);
	 end
	 if ~isempty(outhandles.PA5R)
		disp('...closing PA5R')
		 outhandles.PA5R.status = pa5.closeFunc(PA5R);
	 end
    disp('...closing iodev')
    outhandles.iodev.status = iodev.closeFunc(iodev);
    % Reset TDTINIT
    TDTINIT = 0;
    save(config.TDTLOCKFILE, 'TDTINIT');
    outflag = 1;
else
    disp([mfilename ': TDTINIT is not set!'])
    outflag = -1;
    outhandles = [];
end

%--------------------------------------------------------------------------
% internal function
%--------------------------------------------------------------------------
function [iodev, pa5 ] = setiodev
	iodev.closeFunc = @RPclose; 
	pa5.closeFunc  = @(varargin) -1; 

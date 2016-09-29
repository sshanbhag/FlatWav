function [ outhandles, outflag ] = TDTopen(config, varargin)
%------------------------------------------------------------------------
% [ outhandles, outflag ] = TDTopen(config, varargin)
%------------------------------------------------------------------------
% 
%--- Initializes TDT I/O Hardware ----------------------------------------
% 
%------------------------------------------------------------------------
% Input Arguments:
%   config        matlab struct containing configuration information
%   varargin		for future use
% Output Arguments:
%   outhandles  handle containing indev, outdev, zBUS, PA5L, PA5R
%   outflag    flag to show if TDT is successfully initialized 
%              -1: error 
%               0: not initialized 
%               1: success 
%               2: already initialized
%------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Sharad Shanbhag 
% sshanbhag@neomed.edu
%--------------------------------------------------------------------------
% Created 29 September, 2016 (SJS) from opto_TDTopen.m
%------------------------------------------------------------------------
disp([mfilename ': ...starting TDT devices...']);
outflag = 0; % not initialized

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TDTINIT_FORCE is usually 0, unless user chooses 'RESTART' 
% if TDTINIT is set in the .tdtlock.mat file 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
TDTINIT_FORCE = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% check if the TDT lock file (.tdtlock.mat) exists 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist(config.TDTLOCKFILE, 'file')
    disp([mfilename ': TDT lock file not found: ', config.TDTLOCKFILE]);
    disp('Creating lock file, assuming TDT hardware is not initialized');
    TDTINIT = 0;
    save(config.TDTLOCKFILE, 'TDTINIT');
else
    load(config.TDTLOCKFILE);     % load the lock information
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% check the lock variable (TDTINIT) in the TDT lock file 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if TDTINIT
    questStr = {'TDTINIT already set in .tdtlock.mat', ...
                'TDT Hardware might be active.', ...
                'Continue, Restart TDT Hardware, or Abort?'}; %#ok<UNRCH>
    titleStr = 'FOCHS: TDTINIT error';
    respStr = questdlg(questStr, titleStr, 'Do_Nothing', 'Restart', 'Abort', 'Abort');
    
    switch upper(respStr)
        case 'DO_NOTHING'
            disp([mfilename ': continuing anyway...'])
            outhandles = [];
            outflag = 2;  % already initialized
            return;        
        case 'ABORT'
            disp([mfilename ': aborting initialization...'])
            outhandles = [];
            outflag = -1;  % error state
            return;
        case 'RESTART'
            disp([mfilename ': forcing to start TDT hardware...'])
            TDTINIT_FORCE = 1;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if TDTINIT is not set (TDT hardware not initialized) OR if
% TDTINIT_FORCE is set, initialize TDT hardware
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~TDTINIT || TDTINIT_FORCE
    disp([mfilename ': Configuration = ' config.CONFIGNAME]); 

    %------------------------------------------------------------------
    % initialize the outhandles structure
    %------------------------------------------------------------------
    outhandles = struct();
	 if nargin == 1
	    outhandles.iodev = config.iodev; % Fs, Dnum, Circuit_Path, Circuit_Name
	 else
		 outhandles.iodev = varargin{1};
	 end
    outhandles.iodev.C = [];
    outhandles.iodev.handle = [];
    outhandles.iodev.status = 0;

    %------------------------------------------------------------------
    % setting iodev/PA5 structure (containig function handles)
    % accoring to the config info. Note: setiodev() is defined at the end 
    %------------------------------------------------------------------
    [iodev, pa5] = setiodev;
    
	%------------------------------------------------------------------
	% initialize RX*, PA5; then load and start circuits
	%------------------------------------------------------------------
	try
		% Initialize RX* for input/output
		disp(['...starting ' outhandles.iodev.hardware ' for spike input...'])
		tmpdev = iodev.initFunc('GB', outhandles.iodev.Dnum);
		outhandles.iodev.C = tmpdev.C;
		outhandles.iodev.handle = tmpdev.handle;
		outhandles.iodev.status = tmpdev.status;
		% Initialize Attenuators
		if isempty(pa5)
			disp('using RZ6 internal atten...')
			outhandles.PA5L = [];
			outhandles.PA5R = [];
		else
			disp('...starting PA5...')
			outhandles.PA5L = pa5.initFunc('GB', 1);
			outhandles.PA5R = pa5.initFunc('GB', 2);
		end
		% Load circuits
		disp('...loading circuits...')
		outhandles.iodev.status  = iodev.loadFunc(outhandles.iodev);
	% Starts Circuits
		disp('...starting circuits...')
		iodev.runFunc(outhandles.iodev);
		% Get the input and output sampling rates
		outhandles.iodev.Fs  = iodev.samplefreqFunc(outhandles.iodev);
		disp(['iodev frequency (Hz) = '  num2str(outhandles.iodev.Fs)]);
		% Set the lock
		TDTINIT = 1;  %#ok<NASGU>
		outflag = 1; % success
	catch ME
		TDTINIT = 0; %#ok<NASGU>
		outflag = -1; %#ok<NASGU> % TDT initialization failed  
		disp([mfilename ': error starting TDT hardware'])
		disp(ME.identifier);
		rethrow(ME);
	end

    % save TDTINIT in lock file
    save(config.TDTLOCKFILE, 'TDTINIT');
    return;
end

%--------------------------------------------------------------------------
% internal function
%--------------------------------------------------------------------------
function [iodev, pa5] = setiodev
	iodev.initFunc = @RZ6init; 
	iodev.loadFunc = @RPload;
	iodev.runFunc  = @RPrun;
	iodev.samplefreqFunc = @RPsamplefreq;
	pa5 = [];			

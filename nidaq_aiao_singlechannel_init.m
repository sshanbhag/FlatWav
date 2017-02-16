function NI = nidaq_aiao_singlechannel_init(iface, Dnum)
%--------------------------------------------------------------------------
% NI = nidaq_aiao_singlechannel_init.m
%--------------------------------------------------------------------------
% TytoLogy -> Calibration -> FlatWav program
%--------------------------------------------------------------------------
% initializes nidaq system
%------------------------------------------------------------------------
% Input Arguments:
% 	iface		must be:
% 				'NI' (for traditional interface) 
% 					-OR-
% 				'NI-SESSION' (for new, session interface)
%	Dnum		device id (usually 'Dev1')
% 
% Output Arguments:
% 	NI		struct containing settings for requested type
%		SESSION:
%			NI.S		DAQ Toolbox Session interface object
%	 		NI.chO	analog output channel object
%	 		NI.chI	analog input channel object
%		LEGACY:
%			NI.ao		analog output object (LEGACY only)
%			NI.ai		analog input object (LEGACY only)
%			NI.chO	analog output channel object
%			NI.chI	analog input channel object
%------------------------------------------------------------------------
% See also: NICal, nidaq_aiao_singlechannel_init, DAQ Toolbox (Matlab)
%------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Sharad J Shanbhag
% sshanbhag@neomed.edu
%--------------------------------------------------------------------------
% Created: 4 Dec 2012 (SJS)
% 				Created from nidaq_aiao_init.m
% 
% Revisions:
%	16 Feb 2017 (SJS): updating for Session NI-DAQ interface
%--------------------------------------------------------------------------

%------------------------------------------------------------------------
%------------------------------------------------------------------------
% some checks
%------------------------------------------------------------------------
%------------------------------------------------------------------------
disp('...starting NI hardware...');
if ~any(strcmpi(iface, {'NI', 'NI-SESSION'}))
	errordlg(sprintf('%s: invalid interface %s', mfilename, iface));
	error('%s: invalid interface %s', mfilename, iface);
end

%------------------------------------------------------------------------
%------------------------------------------------------------------------
% Now, Initialize the NI board (PCIe-6351)
%------------------------------------------------------------------------
%------------------------------------------------------------------------

%------------------------------------------------------------------------
%------------------------------------------------------------------------
% Legacy Interface
%------------------------------------------------------------------------
%------------------------------------------------------------------------
if strcmpi(iface, 'NI')
	%---------------------------------------------------------------------
	% 'nidaq' specifies the national instruments device with traditional
	% DAQ Toolbox interface, Device number 1 (get this from the 
	% NI Measurement & Automation Explorer (a.k.a., MAX) program)
	%---------------------------------------------------------------------
	%------------------------------------------------------------------------
	% CONFIGURE ANALOG INPUT SUBSYSTEM
	%------------------------------------------------------------------------
	fprintf('Initializing NIDAQ device for analog input...')
	try
		ai = analoginput('nidaq', Dnum);
		fprintf('...done\n')
	catch errEvent
		fprintf('\nProblem while initializing NIDAQ device!\n\n')
		disp(errEvent)
		return
	end
	% create AI channel
	fprintf('creating analog input channel...')
	chI = addchannel(ai, 0);
	ai.Channel(1).ChannelName = 'responseL';
	fprintf('...done\n');
	%------------------------------------------------------------------------
	% CONFIGURE ANALOG OUTPUT SUBSYSTEM
	%------------------------------------------------------------------------
	fprintf('Initializing NIDAQ device for analog output...')
	try
		ao = analogoutput('nidaq', 'Dev1');
		fprintf('...done\n')
	catch errEvent
		fprintf('\nProblem while initializing NIDAQ device for output!\n\n')
		disp(errEvent)
		return
	end
	% create AO channel
	fprintf('creating analog output channel...')
	chO = addchannel(ao, 0);
	ao.Channel(1).ChannelName = 'stimulusL';
	fprintf('...done\n');
	%-------------------------------------------------------
	% set logging mode
	%	'Disk'	sets logging mode to a file on disk (specified by 'LogFileName)
	%	'Memory'	sets logging mode to memory only
	%	'Disk&Memory'	logs to file and memory
	%-------------------------------------------------------
	set(ai, 'LoggingMode', 'Memory');
	%-------------------------------------------------------
	% save in NI struct
	%-------------------------------------------------------
	NI.ai = ai;
	NI.ao = ao;
	NI.chI = chI;
	NI.chO = chO;
	NI.status = 1;
	return

%------------------------------------------------------------------------
%------------------------------------------------------------------------
% Session Interface
%------------------------------------------------------------------------
%------------------------------------------------------------------------
elseif strcmpi(iface, 'NI-SESSION')
	%---------------------------------------------------------------------
	% CREATE SESSSION
	%---------------------------------------------------------------------
	% 'ni' specifies the national instruments device session
	% Device ID can be found using the daq.getDevices command
	%---------------------------------------------------------------------
	fprintf('Creating session interface...')
	try
		NI.S = daq.createSession('NI');
		fprintf('...done\n');
	catch errEvent
		fprintf('\nCould not create session\n\n');
		disp(errEvent);
		error('Could not create event');
	end
	%---------------------------------------------------------------------
	% CONFIGURE ANALOG INPUT SUBSYSTEM
	%---------------------------------------------------------------------
	fprintf('Adding analog input channels...')
	try
		NI.chI(1) = addAnalogInputChannel(NI.S, Dnum, 0, 'Voltage');
		fprintf('...done\n');
	catch errEvent
		fprintf('\nProblem while adding ai channels to NIDAQ device!\n\n')
		disp(errEvent)
		return
	end
	%---------------------------------------------------------------------
	% CONFIGURE ANALOG OUTPUT SUBSYSTEM
	%---------------------------------------------------------------------
	fprintf('Adding analog output channels...')
	try
		NI.chO(1) = addAnalogOutputChannel(NI.S, Dnum, 0, 'Voltage');
		fprintf('...done\n');
	catch errEvent
		fprintf('\nProblem while adding ao channels to NIDAQ devicet!\n\n')
		disp(errEvent)
		return
	end
	% set outputs to 0
	outputSingleScan(NI.S, zeros(1, length(NI.chO)));

%------------------------------------------------------------------------
%------------------------------------------------------------------------
% Unknown, so throw error
%------------------------------------------------------------------------
%------------------------------------------------------------------------
else
	errordlg(sprintf('Unknown interface %s', iface));
	error('%s: Unknown interface %s', mfilename, iface);
end

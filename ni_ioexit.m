function iodev = ni_ioexit(iodev, debug)
%--------------------------------------------------------------------------
% iodev = ni_ioexit(iodev)
%--------------------------------------------------------------------------
%
% closes NI devices nicely
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Sharad J Shanbhag
% sshanbhag@neomed.edu
%--------------------------------------------------------------------------
% Created: 15 October 2012 (SJS)
% 				Created from NICal_NIexit.m
% 
% Revisions:
%--------------------------------------------------------------------------



%-----------------------------------------------------------------------
%-----------------------------------------------------------------------
% Clean up the RP circuits
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

% save settings information to mat file
if exist('debug', 'var')
	save(fullfile(pwd, 'NI_EventLogs.mat'), ...
			'EventLogAI'			, ...
			'EventLogAO'			, ...
			'-MAT' );
end


	

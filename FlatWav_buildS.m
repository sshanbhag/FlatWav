function S = FlatWav_buildS
%------------------------------------------------------------------------
% S = FlatWav_buildS
%------------------------------------------------------------------------
% 
% Builds synthesis parameters interface struct S.
% 
%------------------------------------------------------------------------
% Input Arguments:
%	none
% 
% Output Arguments:
% 	S
%
%------------------------------------------------------------------------
% See also: FlatWav
%------------------------------------------------------------------------

%------------------------------------------------------------------------
% Sharad J. Shanbhag
% sshanbhag@neomed.edu
%------------------------------------------------------------------------
% Created: 10 Oct, 2011 (SJS)
%
% Revisions:
%------------------------------------------------------------------------
% TO DO:
%------------------------------------------------------------------------

% define synth things

% Sampling rate (default)
S.Fs = 44100;

% Max possible # of synt parameters on GUI
S.MaxShParam = 6;
% generate tag names for the GUI elements.  
%	Text tags are of format 'p[1-6]Text'
%	Ctrl tags are of format 'p[1-6]Ctrl'
for n = 1:S.MaxShParam
	S.CtrlTags{n} = sprintf('p%dCtrl', n);
	S.TextTags{n} = sprintf('p%dText', n);
end

% # of signal types
S.Ntypes = 3;
% names of signal types
S.Types = {'tone', 'noise', 'sweep'};
% Gui Text label values
S.Text =	{	{'Dur', 'Freq', 'Amp', 'RampTime'}; ...
							{'Dur', 'Fmin', 'Fmax', 'Amp', 'RampTime'}; ...
							{'Dur', 'Fmin', 'Fmax', 'Amp', 'RampTime'}; ...
					};
% synth struct parameter names associated with labels/controls
% (stimulus parameters)
S.Param =	{	{'dur', 'freq', 'amp', 'ramp'}; ...
							{'dur', 'fmin', 'fmax', 'amp', 'ramp'}; ...
							{'dur', 'fmin', 'fmax', 'amp', 'ramp'}; ...
					};
% default values
S.DefaultVals = ...
					{	[100, 1000, 1, 5]; ...
						[100, 500, 5000, 1, 5]; ...
						[100, 500, 5000, 1, 5]; ...
					};
% setup list of # of parameters per type
S.Nparam = zeros(S.Ntypes, 1);
for t = 1:S.Ntypes
	S.Nparam(t) = length(S.Text{t});
end

% setup list of ctrl and text handles that match up with Text/Param cells
for t = 1:S.Ntypes
	for n = 1:S.Nparam(t)
		S.CtrlHandles{t}{n} = sprintf('p%dCtrl', n);
		S.TextHandles{t}{n} = sprintf('p%dText', n);
	end
end
				

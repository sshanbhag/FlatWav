function S = FlatWav_buildS
%------------------------------------------------------------------------
% S = FlatWav_buildS
%------------------------------------------------------------------------
% 
% Builds synthesis parameters interface struct S.
% 
%             Fs: 44100
%       MaxNParam: 6
%        CtrlTags: {'p1Ctrl'  'p2Ctrl'  'p3Ctrl'  'p4Ctrl'  'p5Ctrl'  'p6Ctrl'}
%        TextTags: {'p1Text'  'p2Text'  'p3Text'  'p4Text'  'p5Text'  'p6Text'}
%          Ntypes: 3
%           Types: {'tone'  'noise'  'sweep'}
%            Text: {3x1 cell}
%           Param: {3x1 cell}
%     DefaultVals: {3x1 cell}
%          Nparam: [3x1 double]
%     CtrlHandles: {{1x4 cell}  {1x5 cell}  {1x5 cell}}
%     TextHandles: {{1x4 cell}  {1x5 cell}  {1x5 cell}}
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
%	7 Dec 2012 (SJS): updated documentation
%------------------------------------------------------------------------
% TO DO:
%------------------------------------------------------------------------

% define synth things

% Sampling rate (default)
S.Fs = 44100;

% Max possible # of synt parameters on GUI
S.MaxNParam = 6;
% generate tag names for the GUI elements.  
%	Text tags are of format 'p[1-6]Text'
%	Ctrl tags are of format 'p[1-6]Ctrl'
for n = 1:S.MaxNParam
	S.CtrlTags{n} = sprintf('p%dCtrl', n);
	S.TextTags{n} = sprintf('p%dText', n);
end

% # of signal types
S.Ntypes = 3;
% names of signal types
S.Types = {'tone', 'noise', 'sweep'};
% Gui Text label values
S.Text =	{	{'Dur (ms)', 'Freq (Hz)', 'Amp (Vpk)', 'RampTime (ms)'}; ...
				{'Dur (ms)', 'Fmin (Hz)', 'Fmax (Hz)', 'Amp (Vpk)', 'RampTime (ms)'}; ...
				{'Dur (ms)', 'Fmin (Hz)', 'Fmax (Hz)', 'Amp (Vpk)', 'RampTime (ms)'}; ...
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
				

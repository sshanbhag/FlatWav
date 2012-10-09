% define synth things

Fs = 44100;

MaxSynthParam = 6;
for n = 1:MaxSynthParam
	ctrlTags{n} = sprintf('p%dCtrl', n);
	textTags{n} = sprintf('p%dText', n);
end

synthTypes = {'tone', 'noise', 'sweep'}
synthText =		{	{'Dur', 'Freq', 'Amp', 'RampTime'}; ...
						{'Dur', 'Fmin', 'Fmax', 'Amp', 'RampTime'}; ...
						{'Dur', 'Fmin', 'Fmax', 'Amp', 'RampTime'}; ...
					};
synthParam =	{	{'dur', 'freq', 'amp', 'ramp'}; ...
						{'dur', 'fmin', 'fmax', 'amp', 'ramp'}; ...
						{'dur', 'fmin', 'fmax', 'amp', 'ramp'}; ...
					};

synthDefaultVals = ...
					{	[100, 1000, 1, 5]; ...
						[100, 500, 5000, 1, 5]; ...
						[100, 500, 5000, 1, 5]; ...
					};

for t = 1:length(synthTypes)
	for n = 1:length(synthText{t})
		synthCtrlHandles{t}{n} = sprintf('p%dCtrl', n);
		synthTextHandles{t}{n} = sprintf('p%dText', n);
	end
end
				
% initialize structs

% type = 'TONE'
% 
% typenum = find(strcmpi(type, synthTypes))
% 

typenum = 1;
tonedefault.type = synthTypes{typenum};
for n = 1:length(synthText{typenum})
	tonedefault.(synthParam{typenum}{n}) = synthDefaultVals{typenum}(n);
end

typenum = 2;
noisedefault.type = synthTypes{typenum};
for n = 1:length(synthText{typenum})
	noisedefault.(synthParam{typenum}{n}) = synthDefaultVals{typenum}(n);
end

typenum = 3;
sweepdefault.type = synthTypes{typenum};
for n = 1:length(synthText{typenum})
	sweepdefault.(synthParam{typenum}{n}) = synthDefaultVals{typenum}(n);
end

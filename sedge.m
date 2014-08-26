if SMOOTHEDGES
	if LOWCUT
		smoothedges_win = SMOOTHEDGES * [LOWCUT corr_frange];
	else
		smoothedges_win = SMOOTHEDGES * corr_frange;
	end	
end


% get midpoints for smooth windows
midpoints = floor(smoothedges_win ./ 2)

% indices for center of smoothing will be given by valid_indices.
% use this to determine indices of Sadj to be smoothed
% check if lowcut?
if LOWCUT
	lcindx = max(lowcutindices);
	sindx{1} = (lcindx - midpoints(1)):(lcindx + midpoints(1));
	sindx{2} = (valid_indices(1) - midpoints(1)):(valid_indices(1) + midpoints(1));
	sindx{3} = (valid_indices(end) - midpoints(2)):(valid_indices(end) + midpoints(2));
else
	sindx{1} = (valid_indices(1) - midpoints(1)):(valid_indices(1) + midpoints(1));
	sindx{2} = (valid_indices(end) - midpoints(2)):(valid_indices(end) + midpoints(2));
end


figure(1)
plot(SdBadj, 'k')

for n = 1:length(sindx)
	spiece{n} = SdBadj(sindx{n});
	hold on
		plot(sindx{n}, spiece{n}, 'r.')
	hold off
end

for n = 1:length(sindx)
	spiece{n} = moving_average(spiece{n}, 7);
	hold on
		plot(sindx{n}, spiece{n}, 'g-');
	hold off
end


for n = 1:length(sindx)
	SdBadj(sindx{n}) = spiece{n};
end


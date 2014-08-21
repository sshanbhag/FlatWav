IOfigure = figure;
set(IOfigure, 'Position', [10 235 973 500]);
% create subplot axes
% 	subplot('Position',[left bottom width height]) creates an axes at the
% 	position specified by a four-element vector. left, bottom, width, and
% 	height are in normalized coordinates in the range from 0.0 to 1.0
pw = 0.19;
ph = 0.375;
x1 = 0.05;
x2 = 0.3;
x3 = 0.55;
x4 = 0.8;
y1 = 0.55;
y2 = 0.09;

rsig =	subplot('Position', [x1	y1	pw	ph]);
rmag =	subplot('Position', [x2	y1	pw	ph]);
rphi =	subplot('Position', [x3	y1	pw	ph]);
rspec =	subplot('Position', [x4	y1	pw	ph]);

asig =	subplot('Position', [x1	y2	pw	ph]);
amag =	subplot('Position', [x2	y2	pw	ph]);
aphi =	subplot('Position', [x3	y2	pw	ph]);
aspec =	subplot('Position', [x4	y2	pw	ph]);



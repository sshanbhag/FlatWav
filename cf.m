winwidth = 250;
winheight = 150;

%% setup
% open figure
F = figure;
% hide menubar
set(F, 'MenuBar', 'none')
% get position, update width
pos = get(F, 'Position');
pos = [pos(1:2) winwidth winheight];
set(F, 'Position', pos);
% title


% Create the button group.
h = uibuttongroup('visible','off', 'Units', 'points', 'Position', [5 5 175 100]);
% Create three radio buttons in the button group.
u0 = uicontrol('Style','Radio','String','Sound Card (winsound)',...
    'pos',[25 60 150 35],'parent',h,'HandleVisibility','off');
u1 = uicontrol('Style','Radio','String','NI Board (NIDAQ)',...
    'pos',[25 20 150 35],'parent',h,'HandleVisibility','off');
% Initialize some button group properties. 
% set(h,'SelectionChangeFcn',@selcbk);
set(h,'SelectedObject',[]);  % No selection
set(h, 'Title', 'Select Output Device:')
set(h, 'FontSize', 10);
set(h,'Visible','on');


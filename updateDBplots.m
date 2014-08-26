function updateDBplots(hObject, eventdata, handles)


	if any(	[	isempty(handles.respFs) ...
					(isempty(handles.rawresp) && isempty(handles.adjresp))	])
		error('No data for db analysis');
	end
	
	% if dB figure has not been created, do so
	if isempty(handles.dBFigure) || ~ishandle(handles.dBFigure) ...
				|| ~ishandle(handles.P.rawdb) || ~ ishandle(handles.P.adjdb)
		fprintf('Creating dBFigure\n\n');
		handles.dBFigure = figure;
	% 	set(handles.dBFigure, 'Position', [10 235 973 500]);
		% create subplot axes
		handles.P.rawdb =	subplot(211);
		handles.P.adjdb = subplot(212);
		guidata(hObject, handles);
	end

	if ~isempty(handles.respFs)
		rmsbins = ms2samples(handles.PeakRMSWindow, handles.respFs);
	else
		error('updateDBplots: handles.respFs is empty!\n\n');
	end

	% compute rms of raw response, plot
	if ~isempty(handles.rawresp)
		[rawrms, startbins, endbins] = block_rms(handles.rawresp, rmsbins);
		% find peak and peak index of rms values
		[handles.rawrespmax.val, handles.rawrespmax.indx] = max(rawrms);
		% compute peak dB SPL
		rawdBSPL = dbspl(handles.VtoPa*handles.rawrespmax.val);
		% update display
		dbtext = sprintf('Peak Raw dB SPL: %.2f\n', rawdBSPL);
		fprintf('%s\n', dbtext);
		
		% find max point (in milliseconds)
		xval = rmsbins * handles.rawrespmax.indx + (rmsbins ./ 2);
		xval = fix(bin2ms(xval, handles.respFs));
		
		% plot
		subplot(handles.P.rawdb)

		tvec = bin2ms( (1:length(handles.rawresp))-1, handles.respFs);
		
		nrms = length(rawrms);
		x = zeros(2*nrms, 1);
		y = zeros(2*nrms, 1);
		for n = 1:nrms
			x(2 * (n - 1) + 1) = startbins(n);
			x(2 * (n - 1) + 2) = endbins(n);
			y(2 * (n - 1) + 1) = rawrms(n);
			y(2 * (n - 1) + 2) = rawrms(n);
		end
		x = bin2ms(x, handles.respFs);
		y = dbspl(handles.VtoPa*y);
		yresp = max(y) * normalize(handles.rawresp);
		plot(handles.P.rawdb, tvec, yresp, 'b-');
		hold(handles.P.rawdb, 'on')
			plot(handles.P.rawdb, x, y, 'k-', 'LineWidth', 2);
			plot(handles.P.rawdb, xval, rawdBSPL, 'b*', 'MarkerSize', 10);
		hold(handles.P.rawdb, 'off')
		text(xval, rawdBSPL, sprintf('  %.2f', rawdBSPL));
		ylabel(handles.P.rawdb, 'Raw', 'Color', 'r')
		set(handles.P.rawdb, 'XTickLabel', []);
		xlim(handles.P.rawdb, [min(tvec) max(tvec)])
		ylim(handles.P.rawdb, [min(yresp) 1.05*max(y)])
		grid(handles.P.rawdb, 'on')
		% get ticks
		time_ticks = get(handles.P.rawdb, 'XTick');	
		guidata(hObject, handles);
	end
	
		% compute rms of adj response, plot
	if ~isempty(handles.adjresp)
		[adjrms, startbins, endbins] = block_rms(handles.adjresp, rmsbins);
		% find peak and peak index of rms values
		[handles.adjrespmax.val, handles.adjrespmax.indx] = max(adjrms);
		% compute peak dB SPL
		adjdBSPL = dbspl(handles.VtoPa*handles.adjrespmax.val);
		% update display
		dbtext = sprintf('Peak Adj dB SPL: %.2f\n', adjdBSPL);
		fprintf('%s\n', dbtext);
		
		% find max point (in milliseconds)
		xval = rmsbins * handles.rawrespmax.indx + (rmsbins ./ 2);
		xval = fix(bin2ms(xval, handles.respFs));
		
		% plot
		subplot(handles.P.adjdb)

		tvec = bin2ms( (1:length(handles.adjresp))-1, handles.respFs);
		
		nrms = length(adjrms);
		x = zeros(2*nrms, 1);
		y = zeros(2*nrms, 1);
		for n = 1:nrms
			x(2 * (n - 1) + 1) = startbins(n);
			x(2 * (n - 1) + 2) = endbins(n);
			y(2 * (n - 1) + 1) = adjrms(n);
			y(2 * (n - 1) + 2) = adjrms(n);
		end
		x = bin2ms(x, handles.respFs);
		y = dbspl(handles.VtoPa*y);
		yresp = max(y) * normalize(handles.adjresp);
		plot(handles.P.adjdb, tvec, yresp, 'r-');
		hold(handles.P.adjdb, 'on')
			plot(handles.P.adjdb, x, y, 'k-', 'LineWidth', 2);
			plot(handles.P.adjdb, xval, adjdBSPL, 'r*', 'MarkerSize', 10);
		hold(handles.P.adjdb, 'off')
		text(xval, adjdBSPL, sprintf('  %.2f', adjdBSPL));
		ylabel(handles.P.adjdb, 'Adj', 'Color', 'r')
		set(handles.P.adjdb, 'XTickLabel', []);
		xlim(handles.P.adjdb, [min(tvec) max(tvec)])
		ylim(handles.P.adjdb, [min(yresp) 1.05*max(y)])
		grid(handles.P.adjdb, 'on')
		% get ticks
		time_ticks = get(handles.P.adjdb, 'XTick');	
		guidata(hObject, handles);
	end
	
% 		rawresp = handles.rawresp;
% 		respFs = handles.respFs;
% 		adjresp = handles.adjresp;
% 		
%  		save rmsvals.mat rawresp adjresp respFs rmsbins -MAT
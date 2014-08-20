function UpdateInputFilter(hObject, eventdata, handles)
%--------------------------------------------------
%--------------------------------------------------
% Filter settings
% Define a bandpass filter for 
% processing the input data
%--------------------------------------------------
% Nyquist frequency
fnyq = handles.S.Fs / 2;
handles.LPFc = fnyq - 100;
handles.HPFc = 200;
handles.FilterOrder = 4;
% passband definition
fband = [handles.HPFc handles.LPFc] ./ fnyq;
% filter coefficients using a butterworth bandpass filter
[handles.fcoeffb, handles.fcoeffa] = ...
									butter(handles.FilterOrder, fband, 'bandpass');
guidata(hObject, handles);



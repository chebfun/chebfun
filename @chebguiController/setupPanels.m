function handles = setupPanels(handles)
%SETUPPANELS    Populate the panels on CHEBGUI

% Create a toggle button for whether the latest solution should be used as an
% initial guess:
handles.toggle_useLatest = uicontrol('Parent', handles.panel_input, ...
    'Style', 'togglebutton', ...
    'String','Use latest solution as initial guess', ...
    'FontSize', 12, ...
    'Callback', @(hObject, eventdata) ...
        toggleUseLatestCallback(hObject, guidata(hObject)), ...
    'Units', 'normalized', 'Position', [0.1 0.025 0.8 0.075]) ;

end

function toggleUseLatestCallback(hObject, handles)
% Called when user toggles between using the latest solution as an initial
% guess.

% Which state is the button in?
newVal = get(hObject, 'Value');

if ( newVal ) % User wants to use latest solution
    set(handles.input_GUESS, 'String', 'Using latest solution');
else
    set(handles.input_GUESS, 'String', '');
    set(handles.input_GUESS, 'Enable', 'On');
    handles.guifile.init = '';
end

guidata(hObject, handles);

end
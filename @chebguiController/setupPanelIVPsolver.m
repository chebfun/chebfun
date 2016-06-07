function handles = setupPanelIVPsolver(handles)
%SETUPPANELIVPSOLVER   Setup panel for the IVP solver option.

% Background colour for text fields:
textBackgroundColour = get(handles.panel_buttons, 'BackgroundColor');
textFontsize = 12;
leftMargin = 0.025;
textHeight = 0.45;
vertLoc = linspace(.55, .05, 2);

%% Create a panel for the IVP solver option:
handles.panel_IVPsolver = uibuttongroup('parent', handles.panel_buttons, ...
    'Title', 'IVP solver', 'BackgroundColor', textBackgroundColour, ...
    'Position', [0.045 .465 .92 .14], 'FontSize', textFontsize, ...
    'BorderType', 'etchedin', 'Visible', 'off', ...
    'SelectionChangeFcn', @(hObject, eventdata) ...
    ivpSolverSelection(hObject, eventdata, guidata(hObject)));

%% Populate the IVP solver option panel with two buttons.
handles.button_timestepping = uicontrol('Parent', handles.panel_IVPsolver, ...
    'Style', 'radiobutton', 'BackgroundColor', textBackgroundColour, ...
    'String','Time-stepping', 'FontSize', textFontsize, ...
    'Units', 'normalized', ...
    'Position', [leftMargin vertLoc(1) 1-leftMargin textHeight]);

handles.button_global = uicontrol('Parent', handles.panel_IVPsolver, ...
    'Style', 'radiobutton', 'BackgroundColor', textBackgroundColour, ...
    'String','Global', 'FontSize', textFontsize, ...
    'Units', 'normalized', ...
    'Position', [leftMargin vertLoc(2) 1-leftMargin textHeight]);

end

function ivpSolverSelection(hObject, eventdata, handles)
% hObject    handle to the selected object in panel_IVPsolver
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was 
%             selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)

% What's the new IVP solver selected?
newIVPsolver = get(eventdata.NewValue, 'String');

% Change the ivpSolver option stored in handles.guifile.options based on the new
% selection.
if ( strcmp(newIVPsolver, get(handles.button_timestepping, 'String')) )
    % User selected time-stepping. Default time stepping method is ODE113:
    handles.guifile.options.ivpSolver = 'ode113';
    % Hide the initial guess panel:
    set(handles.panel_initialGuess, 'Visible', 'off')
else
    % User selected global solver. Default global method is collocation:
    handles.guifile.options.ivpSolver = 'values';
    % Make the initial guess panel visible:
    set(handles.panel_initialGuess, 'Visible', 'on')
end

% Update HOBJECT.
guidata(hObject, handles);

end
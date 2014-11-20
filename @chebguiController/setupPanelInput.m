function handles = setupPanelInput(handles)
%SETUPPANELS    Populate the panels on CHEBGUI

% Background colour for text fields:
textBackgroundColour = get(handles.panel_input, 'BackgroundColor');
textFontsize = 12;
textHeight = 0.05;
panelLeftMargin = 0.01;
panelWidth = .975;
panelHeight = 0.29;
bottomPanelMargin = 0.005;
inputBoxLeftMargin = 0.025;
inputBoxWidth = .95;


% Create strings above input boxes:
handles.text_LBCs = uicontrol('Parent', handles.panel_input, ...
    'Style', 'text', 'BackgroundColor', textBackgroundColour, ...
    'String','Left boundary condition(s)', 'FontSize', textFontsize, ...
    'Units', 'normalized', 'Position', [0 0.52 1 textHeight], ...
    'Visible', 'off');
handles.text_RBCs = uicontrol('Parent', handles.panel_input, ...
    'Style', 'text', 'BackgroundColor', textBackgroundColour, ...
    'String','Right boundary condition(s)', 'FontSize', textFontsize, ...
    'Units', 'normalized', 'Position', [0 0.39 1 textHeight], ...
    'Visible', 'off');

%% Panel for domain
handles.panel_domain = uipanel('Parent', handles.panel_input, ...
    'Title', 'Domain', 'Titleposition', 'centertop',...
    'BackgroundColor', textBackgroundColour, ...
    'Position', [panelLeftMargin 3*panelHeight + bottomPanelMargin .5 .11], 'FontSize', textFontsize, ...
    'BorderType', 'etchedin');

handles.input_domain = uicontrol('Parent', handles.panel_domain, ...
    'Style', 'edit', ...
    'String','[-1 1]', ...
    'FontSize', 10, 'BackgroundColor', [1 1 1], 'FontName', 'Monospaced', ...
    'Callback', @(hObject, eventdata) ...
        input_domain_Callback(hObject, guidata(hObject)), ...
    'Units', 'normalized', 'Position', [.15 0.15 0.7 .8]);

%% Panel for differential equations
handles.panel_DEs = uipanel('Parent', handles.panel_input, ...
    'Title', 'Differential equations', 'Titleposition', 'centertop',...
    'BackgroundColor', textBackgroundColour, ...
    'Position', [panelLeftMargin 2*panelHeight + bottomPanelMargin panelWidth panelHeight],...
    'FontSize', textFontsize, ...
    'BorderType', 'etchedin');

%% Panel for boundary conditions
handles.panel_BCs = uipanel('Parent', handles.panel_input, ...
    'Title', 'Boundary conditions', 'Titleposition', 'centertop',...
    'BackgroundColor', textBackgroundColour, ...
    'Position', [panelLeftMargin panelHeight + bottomPanelMargin panelWidth panelHeight],...
    'FontSize', textFontsize, ...
    'BorderType', 'etchedin');

%% Setup panel for initial guess
handles.panel_initialGuess = uipanel('Parent', handles.panel_input, ...
    'Title', 'Initial guess',  'Titleposition', 'centertop',...
    'BackgroundColor', textBackgroundColour, ...
    'Position', [0.01 bottomPanelMargin panelWidth panelHeight], ...
    'FontSize', textFontsize, ...
    'BorderType', 'etchedin');

handles.text_initial = uicontrol('Parent', handles.panel_initialGuess, ...
    'Style', 'text', 'BackgroundColor', textBackgroundColour, ...
    'String','Initial condition', 'FontSize', textFontsize, ...
    'Units', 'normalized', 'Position', [0 .85 1 .20], 'visible','off') ;

% Input box for initial guess
handles.input_GUESS = uicontrol('Parent', handles.panel_initialGuess, ...
    'Style', 'edit', 'Max', 2, 'Min',0, ...
    'HorizontalAlignment', 'left', ...
    'FontSize', 10, 'BackgroundColor', [1 1 1], 'FontName', 'Monospaced', ...
    'Callback', @(hObject, eventdata) ...
        input_GUESS_Callback(hObject, guidata(hObject)), ...
    'KeyPressFcn',  @(hObject, eventdata) ...
        input_GUESS_KeyPressFcn(hObject, eventdata, handles), ...
    'ButtonDownFcn', @(hObject, eventdata) ...
        ButtonDownFcn(hObject, eventdata, handles), ...
    'Units', 'normalized', 'Position', [inputBoxLeftMargin 0.35 inputBoxWidth .6]);

% Create a toggle button for whether the latest solution should be used as an
% initial guess:
handles.toggle_useLatest = uicontrol('Parent', handles.panel_initialGuess, ...
    'Style', 'togglebutton', ...
    'String','Use latest solution as initial guess', ...
    'FontSize', 12, 'Enable', 'off', ...
    'Callback', @(hObject, eventdata) ...
        toggleUseLatestCallback(hObject, guidata(hObject)), ...
    'Units', 'normalized', 'Position', [0.1 0.05 0.8 0.25]) ;

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

function input_GUESS_Callback(hObject, handles)
% Plot the initial guess/condition when it is entered in the appropriate field.

% Find the string.
newString = cellstr(get(hObject, 'String'));

% Remove tabs
newString = chebguiController.removeTabs(newString);
set(hObject, 'String', newString);

handles.guifile.init = newString;
if ( isempty(newString) || (iscell(newString) && (numel(newString) == 1) && ...
        isempty(newString{1})) )
    handles.init = '';
    axes(handles.fig_sol);
    cla(handles.fig_sol, 'reset');
    guidata(hObject, handles);
    return
end

loadVariables(handles.importedVar)

guidata(hObject, handles);

% Create the independent space variable:
xtTemp = chebfun(@(x) x, str2num(handles.guifile.domain));

% Assign it to the correct variable, either r, x or t.
if ( ~exist('r', 'var') )
    r = xtTemp;
end

if ( ~exist('x', 'var') )
    x = xtTemp;
end

if ( ~exist('t', 'var') )
    t = xtTemp;
end

% Do something more clever with multiline input
str = cellstr(get(hObject, 'String'));
init = [];
for k = 1:numel(str)
    strk = str{k};
    equalSigns = find(strk == '=');
    if ( numel(equalSigns) > 1 )
        error('CHEBFUN:chebguiWindow:initInput', ...
            'Too many equals signs in input.');
    elseif ( numel(equalSigns) == 1 )
        strk = strk(equalSigns+1:end);
    elseif ( numel(str) > 1 )
        error('CHEBFUN:chebguiWindow:initInput', ...
            ['Error constructing initial guess. Input must include the ' ...
             'names of the dependent variables, i.e. be on the form ' ...
             '"u = %s", ...'], strk)
    end

    strk = deblank(vectorize(strk));
    try
        if ( ~isempty(strk) )
            init = [init eval(strk)]; %#ok<AGROW>
        end
    catch ME
        error('CHEBFUN:chebguiWindow:initInput', ME.message)
    end
end

% Plot the initial guess/solution:
handles.init = init;
axes(handles.fig_sol);
plot(handles.init, 'linewidth', 2)
if ( ~isempty(handles.guifile.options.fixYaxisLower) )
    ylim([str2num(handles.guifile.options.fixYaxisLower) ...
        str2num(handles.guifile.options.fixYaxisUpper)]);
end

% Show grid?
if ( handles.guifile.options.grid )
    grid on
end

% Update the figure and the handles.
guidata(hObject, handles);

end

function loadVariables(importedVar)
% Load variables from the workspace to the workspace of the GUI
fNames = fieldnames(importedVar);
for i = 1:length(fNames)
    assignin('caller', fNames{i}, importedVar.(fNames{i}))
end
end

function input_GUESS_KeyPressFcn(hObject, eventdata, handles)
if ( strcmp(eventdata.Key, 'tab') )
    if ( strcmp(eventdata.Modifier, 'shift') )
        uicontrol(handles.input_BC); 
    else
        uicontrol(handles.button_solve);
        set(handles.button_solve, 'selected', 'on');
    end
end
end

function input_GUESS_ButtonDownFcn(hObject, eventdata, handles)

chebguiEdit('chebguiWindow', handles.chebguimainwindow, 'input_GUESS');
input_GUESS_Callback(hObject, eventdata, handles);

end

function input_domain_Callback(hObject, handles)

in = get(hObject, 'String');
input = str2num(in);

% Checks to see if input is not numeric or empty. If so, default left end
% of the domain is taken to be -1.
if ( input(1) >= input(end) )
    warndlg('Empty domain. Default value [-1, 1] used.')
    in = '[-1, 1]';
    set(hObject, 'String', in);
elseif ( isempty(input) || any(isnan(input)) || (length(input) < 2) )
    warndlg('Domain unrecognized. Default value [-1, 1] used.')
    in = '[-1, 1]';
    set(hObject, 'String', in);
elseif ( ~any(strfind(in, '[')) )
    in = ['[' in ']'];
    set(hObject, 'String', in);
end

set(handles.input_GUESS, 'Enable', 'on');
set(handles.toggle_useLatest, 'Value', 0);
set(handles.toggle_useLatest, 'Enable', 'off');

handles.guifile.domain = in;
guidata(hObject, handles);

end
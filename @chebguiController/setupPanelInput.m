function handles = setupPanelInput(handles)
%SETUPPANELS   Populate the input panel on CHEBGUI

% Background colour for text fields. Ensure that we're using the same one as
% found on the button panel.
textBackgroundColour = get(handles.panel_buttons, 'BackgroundColor');
% Specify the font size of inputs with monospaced font.
monospacedFontsize = 12;
% Specify the font size of panels and text with non-monospaced font.
textFontsize = 12;
% Specify the left margin of the sub-panels.
panelLeftMargin = 0.01;
% Specify the width of the sub-panels (such as inputs for differential equation
% and boundary conditions).
panelWidth = .975;
% Specify the height of the sub-panels.
panelHeight = 0.29;
% The margin where we start the bottom most panel.
bottomPanelMargin = 0.005;
% The left margin for the input boxes inside the subpanels.
inputBoxLeftMargin = 0.025;
% The width of the input boxes within the subpanels.
inputBoxWidth = .95;
%% Setup the input panel
handles.panel_input = uipanel('Parent', handles.mainWindow, ...
    'BackgroundColor', textBackgroundColour, ...
    'Position', [.46305931321540095 .244 .36 .75275], ...
    'BorderType', 'etchedin');

%% Panel for space domain
handles.panel_domain = uipanel('Parent', handles.panel_input, ...
    'Title', 'Domain', 'Titleposition', 'centertop',...
    'BackgroundColor', textBackgroundColour, ...
    'Position', [panelLeftMargin 3*panelHeight+bottomPanelMargin .45 .11], ...
    'FontSize', textFontsize, ...
    'BorderType', 'etchedin');

handles.input_domain = uicontrol('Parent', handles.panel_domain, ...
    'Style', 'edit', 'String','', ...
    'FontSize', monospacedFontsize, 'BackgroundColor', [1 1 1], ...
    'FontName', 'Monospaced', ...
    'Callback', @(hObject, eventdata) ...
        input_domain_Callback(hObject, guidata(hObject)), ...
    'Units', 'normalized', 'Position', [.15 0.15 0.7 .8]);

%% Panel for time interval
handles.panel_timedomain = uipanel('Parent', handles.panel_input, ...
    'Title', 'Time interval', 'Titleposition', 'centertop',...
    'BackgroundColor', textBackgroundColour, ...
    'Position', [.525+panelLeftMargin 3*panelHeight+bottomPanelMargin .45 .11], ...
    'FontSize', textFontsize, ...
    'BorderType', 'etchedin','visible','off');

handles.input_timedomain = uicontrol('Parent', handles.panel_timedomain, ...
    'Style', 'edit', 'String','', 'FontName', 'Monospaced', ...
    'FontSize', monospacedFontsize, 'BackgroundColor', [1 1 1], ...
    'Callback', @(hObject, eventdata) ...
        input_timedomain_Callback(hObject, guidata(hObject)), ...
    'Units', 'normalized', 'Position', [.15 0.15 0.7 .8]);

%% Panel for differential equations
handles.panel_DEs = uipanel('Parent', handles.panel_input, ...
    'Title', 'Differential equations', 'Titleposition', 'centertop',...
    'BackgroundColor', textBackgroundColour, ...
    'Position', [panelLeftMargin 2*panelHeight+bottomPanelMargin panelWidth panelHeight],...
    'FontSize', textFontsize, ...
    'BorderType', 'etchedin');

handles.input_DE = uicontrol('Parent', handles.panel_DEs, ...
    'Style', 'edit', 'Max', 2, 'Min',0, ...
    'HorizontalAlignment', 'left', ...
    'FontSize', monospacedFontsize, 'BackgroundColor', [1 1 1], ...
    'FontName', 'Monospaced', ...
    'Callback', @(hObject, eventdata) ...
        input_DE_Callback(hObject, guidata(hObject)), ...
    'KeyPressFcn',  @(hObject, eventdata) ...
        input_DE_KeyPressFcn(hObject, eventdata, guidata(hObject)), ...
    'ButtonDownFcn', @(hObject, eventdata) ...
        input_DE_ButtonDownFcn(hObject, guidata(hObject)), ...
    'Units', 'normalized', 'Position', [inputBoxLeftMargin 0.1 inputBoxWidth .85]);

%% Panel for boundary conditions
handles.panel_BCs = uipanel('Parent', handles.panel_input, ...
    'Title', 'Boundary conditions', 'Titleposition', 'centertop',...
    'BackgroundColor', textBackgroundColour, ...
    'Position', [panelLeftMargin panelHeight+bottomPanelMargin panelWidth panelHeight],...
    'FontSize', textFontsize, ...
    'BorderType', 'etchedin');

handles.input_BC = uicontrol('Parent', handles.panel_BCs, ...
    'Style', 'edit', 'Max', 2, 'Min', 0, ...
    'HorizontalAlignment', 'left', ...
    'FontSize', monospacedFontsize, 'BackgroundColor', [1 1 1], ...
    'FontName', 'Monospaced', ...
    'Callback', @(hObject, eventdata) ...
        input_BC_Callback(hObject, guidata(hObject)), ...
    'KeyPressFcn',  @(hObject, eventdata) ...
        input_BC_KeyPressFcn(hObject, eventdata, guidata(hObject)), ...
    'ButtonDownFcn', @(hObject, eventdata) ...
        input_BC_ButtonDownFcn(hObject, guidata(hObject)), ...
    'Units', 'normalized', 'Position', [inputBoxLeftMargin 0.1 inputBoxWidth .85],...
     'Interruptible', 'off');

%% Panel for left BCs
handles.panel_leftBCs = uipanel('Parent', handles.panel_input, ...
    'Title', 'Left boundary conditions', 'Titleposition', 'centertop',...
    'BackgroundColor', textBackgroundColour, ...
    'Position', [panelLeftMargin 1.5*panelHeight+bottomPanelMargin panelWidth panelHeight/2],...
    'FontSize', textFontsize, ...
    'BorderType', 'etchedin', 'visible','off');
handles.input_LBC = uicontrol('Parent', handles.panel_leftBCs, ...
    'Style', 'edit', 'Max', 2, 'Min', 0, ...
    'HorizontalAlignment', 'left', ...
    'FontSize', monospacedFontsize, 'BackgroundColor', [1 1 1], ...
    'FontName', 'Monospaced', ...
    'Callback', @(hObject, eventdata) ...
        input_LBC_Callback(hObject, guidata(hObject)), ...
    'KeyPressFcn',  @(hObject, eventdata) ...
        input_LBC_KeyPressFcn(hObject, eventdata, guidata(hObject)), ...
    'ButtonDownFcn', @(hObject, eventdata) ...
        input_LBC_ButtonDownFcn(hObject, guidata(hObject)), ...
    'Units', 'normalized', 'Position', [inputBoxLeftMargin 0.1 inputBoxWidth .85],...
    'Interruptible', 'off');

%% Panel for right BCs
handles.panel_rightBCs = uipanel('Parent', handles.panel_input, ...
    'Title', 'Right boundary conditions', 'Titleposition', 'centertop',...
    'BackgroundColor', textBackgroundColour, ...
    'Position', [panelLeftMargin panelHeight+bottomPanelMargin panelWidth panelHeight/2],...
    'FontSize', textFontsize, ...
    'BorderType', 'etchedin', 'visible','off');
handles.input_RBC = uicontrol('Parent', handles.panel_rightBCs, ...
    'Style', 'edit', 'Max', 2, 'Min', 0, ...
    'HorizontalAlignment', 'left', ...
    'FontSize', monospacedFontsize, 'BackgroundColor', [1 1 1], ...
    'FontName', 'Monospaced', ...
    'Callback', @(hObject, eventdata) ...
        input_RBC_Callback(hObject, guidata(hObject)), ...
    'KeyPressFcn',  @(hObject, eventdata) ...
        input_RBC_KeyPressFcn(hObject, eventdata, guidata(hObject)), ...
    'ButtonDownFcn', @(hObject, eventdata) ...
        input_RBC_ButtonDownFcn(hObject, guidata(hObject)), ...
    'Units', 'normalized', 'Position', [inputBoxLeftMargin 0.1 inputBoxWidth .85],...
    'Interruptible', 'off');

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
    'Style', 'edit', 'Max', 2, 'Min', 0, ...
    'HorizontalAlignment', 'left', ...
    'FontSize', monospacedFontsize, 'BackgroundColor', [1 1 1], ...
    'FontName', 'Monospaced', ...
    'Callback', @(hObject, eventdata) ...
        input_GUESS_Callback(hObject, guidata(hObject)), ...
    'KeyPressFcn',  @(hObject, eventdata) ...
        input_GUESS_KeyPressFcn(hObject, eventdata, guidata(hObject)), ...
    'ButtonDownFcn', @(hObject, eventdata) ...
        input_GUESS_ButtonDownFcn(hObject, guidata(hObject)), ...
    'Units', 'normalized', 'Position', [inputBoxLeftMargin 0.35 inputBoxWidth .6], ...
    'Interruptible', 'off');

% Create a toggle button for whether the latest solution should be used as an
% initial guess:
handles.toggle_useLatest = uicontrol('Parent', handles.panel_initialGuess, ...
    'Style', 'togglebutton', ...
    'String','Use latest solution as initial guess', ...
    'FontSize', 12, 'Enable', 'off', ...
    'Callback', @(hObject, eventdata) ...
        toggleUseLatestCallback(hObject, guidata(hObject)), ...
    'Units', 'normalized', 'Position', [0.1 0.05 0.8 0.25]) ;

%% Setup panel for eigenvalue problem options
handles.panel_eigopts = uipanel('Parent', handles.panel_input, ...
    'Title', 'Eigenvalue problem options',  'Titleposition', 'centertop',...
    'BackgroundColor', textBackgroundColour, ...
    'Position', [0.01 bottomPanelMargin panelWidth panelHeight], ...
    'FontSize', textFontsize, ...
    'BorderType', 'etchedin','visible','off');

handles.text_eigenvalues = uicontrol('Parent', handles.panel_eigopts, ...
    'Style', 'text', 'BackgroundColor', textBackgroundColour, ...
    'String','Number of eigenvalues to seek:', ...
    'FontSize', 12, ...
    'HorizontalAlignment', 'right', ...
    'Units', 'normalized', 'Position', [0 0.65 0.8 0.2]);

handles.edit_eigN = uicontrol('Parent', handles.panel_eigopts, ...
    'Style', 'edit', ...
    'String','', ...
    'FontSize', monospacedFontsize, 'BackgroundColor', [1 1 1], ...
    'FontName', 'Monospaced', ...
    'Callback', @(hObject, eventdata) ...
        edit_eigN_Callback(hObject, guidata(hObject)), ...
    'Units', 'normalized', 'Position', [.81 0.675 0.15 .2], ...
    'Interruptible', 'off');

handles.popupmenu_sigma = uicontrol('Parent', handles.panel_eigopts, ...
    'Style', 'popupmenu', ...
    'String',{'Smoothest eigenmodes', 'Largest magnitude', ...
    'Smallest magnitude', 'Largest real part', 'Smallest real part', ...
    'Largest imag. part', 'Smallest imag. part'}, ...
    'FontSize', 12, 'BackgroundColor', [1 1 1], ...
    'Callback', @(hObject, eventdata) ...
        popupmenu_sigma_Callback(hObject, guidata(hObject)), ...
    'Units', 'normalized', 'Position', [.15 0.15 0.7 .2]);

handles.text_eigenvalues = uicontrol('Parent', handles.panel_eigopts, ...
    'Style', 'text', 'BackgroundColor', textBackgroundColour, ...
    'String','Number of eigenvalues to seek:', ...
    'FontSize', 12, ...
    'HorizontalAlignment', 'right', ...
    'Units', 'normalized', 'Position', [0 0.65 0.8 0.2]);

handles.text_sigma = uicontrol('Parent', handles.panel_eigopts, ...
    'Style', 'text', 'BackgroundColor', textBackgroundColour, ...
    'String','Type of eigenmodes to seek:', ...
    'FontSize', 12, ...
    'HorizontalAlignment', 'center', ...
    'Units', 'normalized', 'Position', [0.1 0.35 0.8 0.2]);
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

% Update HOBJECT with the new information.
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
% This method deals switching focus to the next element when TAB is pressed
% within a CHEBGUI input window.
if ( strcmp(eventdata.Key, 'tab') )
    if ( strcmp(eventdata.Modifier, 'shift') )
        uicontrol(handles.input_BC); 
    else
        uicontrol(handles.button_solve);
        set(handles.button_solve, 'selected', 'on');
    end
end
end

function input_GUESS_ButtonDownFcn(hObject, handles)
% This method deals with opening a special input window if the input box is
% right-clicked.
chebguiEdit('chebguiWindow', handles.chebguimainwindow, 'input_GUESS');
input_GUESS_Callback(hObject, handles);

end

function input_domain_Callback(hObject, handles)
% Called when the spacial domain input has been changed.

% Get the input, and checks it validity.
in = get(hObject, 'String');
input = str2num(in);

% Checks to see if input is not numeric or empty. If so, default left end
% of the domain is taken to be -1.
if ( input(1) >= input(end) )
    warndlg('Empty domain. Default value [-1,1] used.')
    in = '[-1, 1]';
    set(hObject, 'String', in);
elseif ( isempty(input) || any(isnan(input)) || (length(input) < 2) )
    warndlg('Domain unrecognized. Default value [-1,1] used.')
    in = '[-1, 1]';
    set(hObject, 'String', in);
elseif ( ~any(strfind(in, '[')) )
    in = ['[' in ']'];
    set(hObject, 'String', in);
end

% When the domain has been changed, we can't expect to be able to use the
% previously found solution as the initial guess.
set(handles.input_GUESS, 'Enable', 'on');
set(handles.toggle_useLatest, 'Value', 0);
set(handles.toggle_useLatest, 'Enable', 'off');

% Update the domain stored in the guifile.
handles.guifile.domain = in;

% Update HOBJECT with the new information.
guidata(hObject, handles);

end

function input_DE_Callback(hObject, handles)
% Called when the differential equation is entered.

% Obtain the input:
str = cellstr(get(hObject, 'String'));

% Remove tabs:
str = chebguiController.removeTabs(str);

% Update the DE input and store in guifile:
set(handles.input_DE, 'String', str);
handles.guifile.DE = str;

% Auto PDE and EIG detection
for k = 1:numel(str)
    strk = str{k};
    if ( any(strfind(strk, '_')) )
        if ( ~get(handles.button_pde, 'value') )
            handles = chebguiController.switchMode(handles, 'pde');
        end
        break
    elseif ( any(strfind(strk, 'lam') | strfind(strk, 'lambda')) )
        if ( ~get(handles.button_eig, 'value') )
            handles = chebguiController.switchMode(handles, 'eig');
        end
        break
    end
end

% Update HOBJECT with the new information.
guidata(hObject, handles);

end

function input_DE_KeyPressFcn(hObject, eventdata, handles)
% This method deals switching focus to the next element when TAB is pressed
% within a CHEBGUI input window.
if ( strcmp(eventdata.Key, 'tab') )
    if ( strcmp(eventdata.Modifier, 'shift') )
        if ( get(handles.button_pde, 'value') )
            uicontrol(handles.input_timedomain); 
        else
            uicontrol(handles.input_domain); 
        end
    elseif ( get(handles.button_pde, 'value') )
        uicontrol(handles.input_LBC); 
    else
        uicontrol(handles.input_BC); 
    end
end
end

function input_DE_ButtonDownFcn(hObject, handles)
% Open an edit box if the input is right-clicked.
chebguiEdit('chebguiWindow', handles.chebguimainwindow, 'input_DE');
input_DE_Callback(hObject, handles);

end

function input_BC_ButtonDownFcn(hObject, handles)
% Open an edit box if the input is right-clicked.
chebguiEdit('chebguiWindow', handles.chebguimainwindow, 'input_BC');
input_BC_Callback(hObject, handles);

end

function input_BC_Callback(hObject, handles)
% Called after the boundary conditions have been entered.
newString = cellstr(get(hObject, 'String'));
newString = chebguiController.removeTabs(newString); % Remove tabs
set(hObject, 'String', newString);
handles = chebguiController.callbackBCs(handles, newString, 'bc');
handles.guifile.BC = newString;

% Check whether we need to switch from BVP mode to IVP mode or vice versa. This
% can happen if we go from specifying conditions at both ends to only one end
% (or vice versa).
if (get(handles.button_bvp,'value') )
    % Is the problem now a BVP, IVP or FVP?
    isIorF = isIVPorFVP(handles.guifile);
    if ( isIorF )
        handles = chebguiController.switchMode(handles, 'ivp');
        handles.guifile.type = 'ivp';
    end
elseif ( get(handles.button_ivp,'value') )
    % Is the problem now a BVP?
    isIorF = isIVPorFVP(handles.guifile);
    if ( ~isIorF )
        handles = chebguiController.switchMode(handles, 'bvp');
        handles.guifile.type = 'bvp';
    end    
end

% Update HOBJECT with the new information.
guidata(hObject, handles);
end

function input_BC_KeyPressFcn(hObject, eventdata, handles)
% This method deals switching focus to the next element when TAB is pressed
% within a CHEBGUI input window.
if ( strcmp(eventdata.Key, 'tab') )
    if ( strcmp(eventdata.Modifier, 'shift') )
        uicontrol(handles.input_DE); 
    else
        uicontrol(handles.input_GUESS); 
    end
end
end

function input_LBC_Callback(hObject, handles)
% Called after the left boundary condition has been entered (in PDE mode).
newString = cellstr(get(hObject, 'String'));
newString = chebguiController.removeTabs(newString); % Remove tabs
set(hObject, 'String', newString);
handles = chebguiController.callbackBCs(handles, newString, 'lbc');
handles.guifile.LBC = newString;
guidata(hObject, handles);
end


function input_RBC_Callback(hObject, handles)
% Called after the right boundary condition has been entered (in PDE mode).
newString = cellstr(get(hObject, 'String'));
newString = chebguiController.removeTabs(newString); % Remove tabs
set(hObject, 'String', newString);
handles = chebguiController.callbackBCs(handles, newString, 'rbc');
handles.guifile.RBC = newString;
guidata(hObject, handles);
end

function input_LBC_KeyPressFcn(hObject, eventdata, handles)
% This method deals switching focus to the next element when TAB is pressed
% within a CHEBGUI input window.
if ( strcmp(eventdata.Key, 'tab') )
    if ( strcmp(eventdata.Modifier, 'shift') )
        uicontrol(handles.input_DE); 
    else
        uicontrol(handles.input_RBC); 
    end
end
end

function input_RBC_KeyPressFcn(hObject, eventdata, handles)
% This method deals switching focus to the next element when TAB is pressed
% within a CHEBGUI input window.
if ( strcmp(eventdata.Key, 'tab') )
    if ( strcmp(eventdata.Modifier, 'shift') )
        uicontrol(handles.input_LBC); 
    else
        uicontrol(handles.input_GUESS); 
    end
end
end

function input_LBC_ButtonDownFcn(hObject, handles)
% Open an edit box if the input is right-clicked.
chebguiEdit('chebguiWindow', handles.chebguimainwindow, 'input_LBC');
input_LBC_Callback(hObject, handles);

end

function input_RBC_ButtonDownFcn(hObject, handles)
% Open an edit box if the input is right-clicked.
chebguiEdit('chebguiWindow', handles.chebguimainwindow, 'input_RBC');
input_RBC_Callback(hObject, handles);

end

function input_timedomain_Callback(hObject, handles)
% Called after the timedomain has been modified.

% Passing time range for PDEs.
str = get(hObject, 'String');
if ( iscell(str) )
    str = str{:};
end
num = str2num(str);

% We'll input dialogs if we're not happy about the input. The following option
% ensures those will be modal.
options.WindowStyle = 'modal';

% Indicates we had timedomain with negative spacing (0:-.1:1)
while ( isempty(num) )
    str = inputdlg(['Time domain should be a vector of length > 2, with ' ...
        'positive spacing, at which times solution is returned'], ['Set ' ...
        'time domain'],  1, {'0:.1:1'}, options);
    if ( isempty(str) )
        str = '';
        break
    end
    str = str{:};
    num = str2num(str);
end

% The time interval only consisted of two entries.
while ( ~isempty(str) && (numel(num) < 3) )
    h = (num(2) - num(1))/20;
    def = sprintf('%s:%s:%s', num2str(num(1), '%0.0f'), num2str(h, '%0.2g'), ...
        num2str(num(2), '%0.0f'));
    str = inputdlg(['Time domain should be a vector of length > 2 at which ' ...
        'times solution is returned'], 'Set time domain',  1, {def}, options);
    if ( isempty(str) )
        str = '';
        break
    end
    str = str{:};
    num = str2num(str);
end

% Update the GUI with the new accepted input.
set(handles.input_timedomain, 'String', str);
handles.guifile.timedomain = str;

% Update the HOBJECT.
guidata(hObject, handles);

end

function edit_eigN_Callback(hObject, handles)
% Called when the number of eigenvalues we look for changes.
in = get(handles.edit_eigN, 'String');

% Ensure that we got an integer input!
if ( ~isempty(in) && isempty(str2num(in)) )
    errordlg('Invalid input. Number of eigenvalues must be an integer.', ...
        'Chebgui error', 'modal');
    set(handles.edit_eigN, 'String', handles.guifile.options.numeigs);
else
    handles.guifile.options.numeigs = in;
end

% Update HOBJECT with the new info.
guidata(hObject, handles);
end

function popupmenu_sigma_Callback(hObject, handles)
% Called when the type of eigenmodes we seek gets changes. This goes from a
% human-friendly description, such as 'smallest magnitude', to an option
% accepted by EIGS, such as 'sm'.
selected = get(hObject, 'Value');
switch ( selected )
    case 1
        handles.guifile.sigma = '';
    case 2
        handles.guifile.sigma = 'lm';
    case 3
        handles.guifile.sigma = 'sm';
    case 4
        handles.guifile.sigma = 'lr';
    case 5
        handles.guifile.sigma = 'sr';
    case 6
        handles.guifile.sigma = 'li';
    case 7
        handles.guifile.sigma = 'si';
end
guidata(hObject, handles);
end

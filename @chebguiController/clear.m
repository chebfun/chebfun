function handles = clear(handles)
%CLEAR   Reset CHEBGUI, clear all fields and reset options.
%  Calling sequence:
%       HANDLES = CLEAR(HANDLES)
%  where
%       HANDLES:    A Matlab handle of the CHEBGUI figure.
%
%  This method is called when the user presses the 'Clear' button on the CHEBGUI
%  figure.

% Copyright 2017 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% Clear the input fields:
set(handles.input_domain, 'String', '');
set(handles.input_timedomain, 'String', '');
set(handles.input_DE, 'String', '');
set(handles.input_DE, 'String', '');
set(handles.input_LBC, 'String', '');
set(handles.input_RBC, 'String', '');
set(handles.input_BC, 'String', '');
set(handles.input_GUESS, 'String', '');
set(handles.menu_tolerance, 'UserData', '5e-13'); % The default tolerance

% Clear the figures:
chebguiController.initialiseFigures(handles)

% Hide the iteration information:
set(handles.iter_list, 'String', '');
set(handles.iter_text, 'Visible', 'Off');
set(handles.iter_list, 'Visible', 'Off');

% Disable export figures buttons:
set(handles.button_figsol, 'Enable', 'off');
set(handles.button_solve, 'Enable', 'on');
set(handles.button_solve, 'String', 'Solve');

% Enable RHS of BCs again:
set(handles.input_LBC, 'Enable', 'on');
set(handles.input_RBC, 'Enable', 'on');

% Reset eigenvalue options:
set(handles.edit_eigN, 'String', '6');
set(handles.popupmenu_sigma, 'Value', 1)

% We don't have a solution available anymore:
handles.hasSolution = 0;
set(handles.button_exportsoln, 'Value', 0);

end


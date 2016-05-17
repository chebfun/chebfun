function handles = initalizeFields(handles)
%INITIALIZEFIELDS    Initialize fonts when CHEBGUI is started.
%
%  Calling sequence:
%       HANDLES = INITIALIZEFIELDS(HANDLES)
%  where
%       HANDLES:    A Matlab handle of the CHEBGUI figure.
%

% Copyright 2016 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% Set the default font and font size used:
myFont = 'Monospaced';
mySize = 12;

% Obtain a list of all GUI elements.
names = fieldnames(handles);
% Find what elements are input and button elements:
inputLocs = strfind(names, 'input_');
buttonLocs = strfind(names, 'button_');
% Also want the same font for the iter_list information box.
iterListLoc = strfind(names, 'iter_list');
% And what fields are text fields:
textLocs = strfind(names, 'text_');
% Combine all the locations of elements whose font we wish to change
monospacedLocs = strcat(inputLocs, iterListLoc);
variableWidthLocs = strcat(buttonLocs, textLocs);


% Loop through the elements we want to specify the font of.
for fieldCounter = 1:length(inputLocs)
    if ( ~isempty(monospacedLocs{fieldCounter}) )
        % Access the field values dynamically using the .( ) call.
        set(handles.(names{fieldCounter}), 'FontName', myFont);
        set(handles.(names{fieldCounter}), 'FontSize', mySize + 1);
    end
    
    if ( ~isempty(variableWidthLocs{fieldCounter}) )
        set(handles.(names{fieldCounter}), 'FontSize', mySize);
    end
end

% Set the string for popup-menu for the choice of plots:
set(handles.popupmenu_bottomFig,'String', ...
    {'Convergence of Newton iteration', ...
    'Coefficients of the solution'}, ...
    'fontsize', handles.fontsizePanels);

end


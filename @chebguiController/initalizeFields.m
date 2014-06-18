function handles = initalizeFields(handles)
%INITIALIZEFIELDS    Initialize fonts when CHEBGUI is started.
%
%  Calling sequence:
%       HANDLES = INITIALIZEFIELDS(HANDLES)
%  where
%       HANDLES:    A Matlab handle of the CHEBGUI figure.
%

% Copyright 2014 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% Set input boxes to command window font:
try
    % Obtain the font information:
    s = char(com.mathworks.services.FontPrefs.getCodeFont);
    idx = strfind(s, 'name=');
    s1 = s(idx+5:end);
    idx = strfind(s1, ',');
    s1 = s1(1:idx(1)-1);
    myFont = s1;
    idx = strfind(s, 'size=');
    s2 = s(idx+5:end-1);
    mySize = str2double(s2) + 1;
catch
    % If this fails, fall back to this default:
    myFont = 'Monospaced';
    mySize = 15;
end

% Obtain a list of all GUI elements.
names = fieldnames(handles);
% Find what elements are input elements
inputLocs = strfind(names, 'input_');
% Also want the same font for the iter_list information box.
iterListLoc = strfind(names, 'iter_list');
% And what fields are text fields:
textLocs = strfind(names, 'text_');
% Combine all the locations of elements whose font we wish to change
allLocs = strcat(inputLocs, iterListLoc);

% Loop through the elements we want to specify the font of.
for fieldCounter = 1:length(inputLocs)
    if ( ~isempty(allLocs{fieldCounter}) )
        % Access the field values dynamically using the .( ) call.
        set(handles.(names{fieldCounter}), 'FontName', myFont);
        set(handles.(names{fieldCounter}), 'FontSize', mySize);
    end
    
    if ( ~isempty(textLocs{fieldCounter}) )
        set(handles.(names{fieldCounter}), 'FontSize', mySize);
    end
end

% Set the string for popup-menu for the choice of plots:
set(handles.popupmenu_bottomFig,'String', ...
    {'Convergence of Newton iteration', ...
    'Chebyshev coefficients'});

% Set a multiline string for the ultraspherical options
set(handles.button_ultraS, 'String', ...
    '<HTML><BODY>Ultra- <br> spherical</BODY> </HTML>')

end


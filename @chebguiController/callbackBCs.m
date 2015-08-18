function handles = callbackBCs(handles, inputString, type)
%CALLBACKBCS   Govern the behaviour of BC input fields on CHEBGUI.
% Calling sequence:
%   HANDLES = CALLBACKBCS(HANDLES, INPUTSTRING, TYPE)
% where
%   HANDLES:        A MATLAB handle object corresponding to the CHEBGUI figure.
%   INPUTRSTRING:   The input from the user to the BC field.
%   TYPE:           Whether the input was passed to the BC field (used in BVP
%                   and EIG modes) or the LBC and RBC fields (used in PDE mode). 

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% For systems we check one row at a time.
if ( ~iscell(inputString) )
    newString = cell(size(inputString));
    for stringCounter = 1:size(inputString, 1)
        newString{stringCounter,:} = inputString(stringCounter,:);
    end
else
    newString = inputString;
end

% See whether we got passed 'periodic' -- if so, disable the other field.
periodicFlag = false; % Periodic flag
for k = 1:numel(newString)
    if ( ~isempty(strfind(newString{k}, '@')) ...
            || strcmpi(newString{k}, 'dirichlet') ...
            || strcmpi(newString{k}, 'neumann') ...
            || ~isempty(str2num(newString{k})) )
        break
    elseif ( strcmpi(newString{k}, 'periodic') )
        periodicFlag = true;
        break
    end
end

if ( strcmp(type, 'lbc') )
    % Deal with the LBC input.
    if ( periodicFlag )
        set(handles.input_RBC, 'String', 'periodic');
        handles.guifile.RBC = 'periodic';
        set(handles.input_RBC, 'Enable', 'off');
    else
        set(handles.input_RBC, 'Enable', 'on');
        % Clear the information about periodic BCs from the CHEBGUI object.
        % The .LBC field will be loaded in the chebguiwindow.m method.
        if ( strcmp(handles.guifile.LBC, 'periodic') )
            handles.guifile.RBC = '';
            set(handles.input_RBC, 'String', '');
        end
    end
    
elseif ( strcmp(type, 'rbc') )
    % Deal with the RBC input.
    if ( periodicFlag )
        set(handles.input_LBC, 'String', 'periodic');
        handles.guifile.LBC = 'periodic';
        set(handles.input_LBC, 'Enable', 'off');
    else
        set(handles.input_LBC, 'Enable', 'on');
        % Clear the information about periodic BCs from the CHEBGUI object.
        % The .RBC field will be loaded in the chebguiwindow.m method.
        if ( strcmp(handles.guifile.RBC, 'periodic') )
            handles.guifile.LBC = '';
            set(handles.input_LBC, 'String', '');
        end
    end
    
elseif ( strcmp(type, 'bc') )
    
else
    error('CHEBFUN:CHEBGUICONTROLLER:callbackBCs:unknown', ...
        'Unknown BC type %s.', type);
    
end

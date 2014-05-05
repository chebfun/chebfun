function handles = callbackBCs(guifile, handles, inputString, type)

% Copyright 2011 by The University of Oxford and The Chebfun Developers. 
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.

% For systems we check a row at a time

flag = false;  % dirichlet/neumann flag
flag2 = false; % periodic flag

if ~iscell(inputString)
    for stringCounter = 1:size(inputString,1)
        newString{stringCounter,:} = inputString(stringCounter,:);
    end
else
    newString = inputString;
end

for k = 1:numel(newString)
    if ~isempty(strfind(newString{k},'@')) || strcmpi(newString{k},'dirichlet') ...
            || strcmpi(newString{k},'neumann') || ~isempty(str2num(newString{k}))
        flag = true; break
    elseif strcmpi(newString{k},'periodic');
        flag = true; flag2 = true; break
    end
end

if strcmp(type,'lbc')   
    if flag2
        set(handles.input_RBC,'String','periodic');
        handles.guifile.RBC = 'periodic';
        set(handles.input_RBC,'Enable','off');
    else
        set(handles.input_RBC,'Enable','on');
        % Clear the information about periodic BCs from the chebgui object.
        % The .LBC field will be loaded in the chebguiwindow.m method.
        if strcmp(handles.guifile.LBC,'periodic')
            handles.guifile.RBC = '';
            set(handles.input_RBC,'String','');
        end
    end
elseif strcmp(type,'rbc') 
    if flag2
        set(handles.input_LBC,'String','periodic');
        handles.guifile.LBC = 'periodic';
        set(handles.input_LBC,'Enable','off');
    else
        set(handles.input_LBC,'Enable','on');
        % Clear the information about periodic BCs from the chebgui object.
        % The .RBC field will be loaded in the chebguiwindow.m method.
        if strcmp(handles.guifile.RBC,'periodic')
            handles.guifile.LBC = '';
            set(handles.input_LBC,'String','');
        end
    end
elseif strcmp(type,'bc') 
    if flag2
        set(handles.input_LBC,'String','periodic');
        handles.guifile.LBC = 'periodic';
        set(handles.input_LBC,'Enable','off');
        set(handles.input_RBC,'String','periodic');
        handles.guifile.RBC = 'periodic';
        set(handles.input_RBC,'Enable','off');
    else
        if strcmp(get(handles.input_LBC,'String'),'periodic');
            handles.guifile.LBC = '';
            set(handles.input_LBC,'String','');
        end
        set(handles.input_LBC,'Enable','on');
        if strcmp(get(handles.input_RBC,'String'),'periodic');
            handles.guifile.RBC = '';
            set(handles.input_RBC,'String','');
        end
        set(handles.input_RBC,'Enable','on');
    end
end
    
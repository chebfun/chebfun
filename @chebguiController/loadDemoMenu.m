function loadDemoMenu(handles)
%LOADDEMOMENU    Populate the 'Demos' menu on the CHEBGUI figure.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Begin by checking whether we have already loaded the demos
if ( isfield(handles,'demosLoaded') )
    return
end

% Find the folders which demos are stored in. The chebguiDemos folder lives in
% the trunk folder, find the path of the Chebfun trunk.
trunkPath = chebfunroot();

% Append directory information
bvppath = fullfile(trunkPath, 'chebguiDemos', 'bvpdemos');
ivppath = fullfile(trunkPath, 'chebguiDemos', 'ivpdemos');
pdepath = fullfile(trunkPath, 'chebguiDemos', 'pdedemos');
eigpath = fullfile(trunkPath, 'chebguiDemos', 'eigdemos');

% Setup BVPs
D = dir(bvppath);
for demoCounter = 1:length(D)
    demoPath = fullfile(bvppath, D(demoCounter, :).name);
    if ( isempty(strfind(demoPath, '.guifile')) )
        % Only want to load files ending in .guifile
        continue
    end

    % Parse the .GUIFILE that stores the demo:
    [demoName, demoFun, demoType] = parseDemoFile(demoPath, handles);
    
    switch demoType
        % Have three categories of ODE demos. The call to UIMENU() below
        % populuates the menus.
        case 'bvp'
            uimenu('Parent', handles.menu_bvps, ...
                'Label', demoName, ...
                'Separator', 'off', ...
                'HandleVisibility', 'callback',  ...
                'Callback', demoFun);
        case 'system'
            uimenu('Parent', handles.menu_systems, ...
                'Label', demoName, ...
                'Separator', 'off', ...
                'HandleVisibility','callback', ...
                'Callback', demoFun);
    end
end

% Setup IVPs
D = dir(ivppath);
for demoCounter = 1:length(D)
    demoPath = fullfile(ivppath, D(demoCounter, :).name);
    if ( isempty(strfind(demoPath, '.guifile')) )
        % Only want to load files ending in .guifile
        continue
    end

    % Parse the .GUIFILE that stores the demo:
    [demoName, demoFun, demoType] = parseDemoFile(demoPath, handles);
    
    switch demoType
        % Have three categories of ODE demos. The call to UIMENU() below
        % populuates the menus.
        case 'scalar'
            uimenu('Parent', handles.menu_ivps, ...
                'Label', demoName, ...
                'Separator', 'off', ...
                'HandleVisibility', 'callback',  ...
                'Callback', demoFun);
        case 'system'
            uimenu('Parent', handles.menu_IVPsystems, ...
                'Label', demoName, ...
                'Separator', 'off', ...
                'HandleVisibility','callback', ...
                'Callback', demoFun);
    end
end

% Setup PDEs
D = dir(pdepath);
for demoCounter = 1:length(D) % First two entries are . and ..
    demoPath = fullfile(pdepath, D(demoCounter,:).name);

    % Only want to load files ending in .guifile
    if ( isempty(strfind(demoPath, '.guifile')) )
        continue
    end

    % Parse the .GUIFILE that stores the demo:
    [demoName, demoFun, demoType] = parseDemoFile(demoPath, handles);
    
    switch demoType
        case 'scalar'
            uimenu('Parent', handles.menu_pdesingle, ...
                'Label', demoName, ...
                'Separator', 'off', ...
                'HandleVisibility', 'callback', ...
                'Callback', demoFun);
        case 'system'
            uimenu('Parent', handles.menu_pdesystems, ...
                'Label', demoName, ...
                'Separator', 'off', ...
                'HandleVisibility', 'callback', ...
                'Callback', demoFun);
    end
end

% Setup EIGs
D = dir(eigpath);
for demoCounter = 1:length(D) % First two entries are . and ..
    demoPath = fullfile(eigpath,D(demoCounter,:).name);
    if isempty(strfind(demoPath,'.guifile')) 
        continue % Only want to load files ending in .guifile
    end

    % Parse the .GUIFILE that stores the demo:
    [demoName, demoFun, demoType] = parseDemoFile(demoPath, handles);
    
    switch demoType
        case 'scalar'
            uimenu('Parent', handles.menu_eigsscalar, ...
                'Label', demoName, ...
                'Separator', 'off', ...
                'HandleVisibility', 'callback', ...
                'Callback', demoFun);
        case 'system'
            uimenu('Parent', handles.menu_eigssystem, ...
                'Label', demoName, ...
                'Separator', 'off', ...
                'HandleVisibility', 'callback', ...
                'Callback', demoFun);
    end
end

% Notify that we have loaded demos to prevent reloading
set(handles.menu_demos, 'UserData', 1);

end

function [demoName, demoFun, demoType] = parseDemoFile(demoPath, handles)
%PARSEDEMOFILE    Obtain information from the .guifile currently considered.

    % Need to obtain the name and type of the demo as well
    fid = fopen(demoPath);

    % Loop through the problem descriptions for testing at the top of the files.
    while ( true )
        tline = fgetl(fid);
        % Once we don't have a # at the start, we longer are in the problem
        % description part of the file.
        if ( isempty(tline) || tline(1) ~= '#')
            break
        end
    end

    % Now load the problem description and type lines, which will be the next
    % lines up.
    
    % Throw away ' at the ends of the string
    demoName = tline;
    demoName = demoName(2:end-1);
    demoType = fgetl(fid);
    demoType = demoType(2:end-1);

    % Close the file
    fclose(fid);

    % When a user selects a demo in the menu, the callback method of the
    % associated menu item gets called. The lines below create a new function
    % that gets assigned as the callback method of the items -- i.e., when the
    % name of the demo is clicked, demoFun() of the corresponding demo gets
    % called.
    demoFun = @(hObject, eventdata) ...
        hOpenMenuitemCallback(hObject, eventdata, guidata(hObject), demoPath);
end

function hOpenMenuitemCallback(hObject, eventdata, handles, demoPath)
%HOPENMENUITEMCALLBACK  The callback method of the 'Demos' menu items.
%
% When a demo is selected, this method is called. It then populates the CHEBGUI
% figure with the demo.

% Callback function run when the Open menu item is selected
handles.guifile = chebgui.demo2chebgui(demoPath);

% Switch the mode of the GUI according to the type of the problem.
chebguiController.switchMode(handles, handles.guifile.type);

% Populate the CHEBGUI figure.
chebguiController.populate(hObject, handles, handles.guifile);

% We no longer have a solution.
handles.hasSolution = 0;
set(handles.button_exportsoln, 'Enable', 'off');

% Update handle structure
guidata(hObject, handles);

end

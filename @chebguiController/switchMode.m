function handles = switchMode(handles, newMode)
%SWITCHMODE   Go from one mode of CHEBGUI to another.
% Calling sequence:
%   HANDLES = SWITCHMODE(HANDLES, NEWMODE)
% where
%   HANDLES:    A MATLAB handle object for the CHEBGUI figure.
%   NEWMODE:    Mode we want to switch to.

% Copyright 2016 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% Do a lot of disabling/enabling and hiding/showing objects on the CHEBGUI
% figure.
if ( strcmp(newMode, 'bvp') ) % Going into BVP mode
    handles.guifile.type = 'bvp';
    
    set(handles.button_bvp, 'Value', 1)
    set(handles.button_ivp, 'Value', 0)
    set(handles.button_pde, 'Value', 0)
    set(handles.button_eig, 'Value', 0)
    
    set(handles.panel_DEs, 'Title', 'Differential equations')
    set(handles.toggle_useLatest, 'Visible', 'on')
    
    set(handles.panel_timedomain, 'Visible', 'off')
    
    set(handles.button_realplot, 'Visible', 'off')
    set(handles.button_imagplot, 'Visible', 'off')
    set(handles.popupmenu_bottomFig, 'Visible', 'on')
    set(handles.toggle_useLatest, 'Value', 0)
    
    if ( isempty(fieldnames(handles.latest)) )
        set(handles.toggle_useLatest, 'Enable', 'off')
    end
    
    % Change heading for BC/IC input field:
    set(handles.panel_BCs, 'Title', 'Boundary conditions')
    
    % Hide the iter_list box:
    handles = hideIterList(handles);
    
    % Show the initial guess panel, and ensure it has the correct title:
    set(handles.panel_initialGuess, 'Visible', 'on', 'Title', 'Initial guess');
    
    set(handles.panel_eigopts, 'Visible', 'Off')
        
    % Switch from LBC and RBC fields to BC
    handles = toggleBCinput(handles, 'off', 'on');
    
    % Change the list of available options.
    % Enable ODE menu options
    set(handles.menu_odedampednewton, 'Enable', 'On')
    set(handles.menu_odeplotting, 'Enable', 'On')
    % Disable IVP solver option
    set(handles.menu_ivpSolver, 'Enable', 'off');
    % Disable PDE menu options
    set(handles.menu_pdeSolver, 'Enable', 'Off')
    set(handles.menu_pdeplotting, 'Enable', 'Off')
    set(handles.menu_pdeholdplot, 'Enable', 'Off')
    set(handles.menu_pdefix, 'Enable', 'Off')
    set(handles.menu_fixN, 'Enable', 'Off')

    % Make the discretization visible
    set(handles.panel_discretization, 'Visible', 'on');
    
    % Hide the IVP solver panel visible
    set(handles.panel_IVPsolver, 'Visible', 'off')
    
    % Clear the top figure
    chebguiController.initialiseFigureTop(handles)

    % Always clear the bottom figure
    chebguiController.initialiseFigureBottom(handles)
        
    % Make the discretization panel visible
    set(handles.panel_discretization, 'Visible', 'on')
    
    % Set correct header to top figure
    set(handles.panel_figSol, 'Title', 'Solution')
    
    % Note: The line below is a nice way to disable to objects, rather than hide
    % them.
    % Enable all childrens of discretisation panel
    %     set(findall(handles.panel_discretization,  '-property',  'enable'),  ...
    %         'visible',  'on')
elseif ( strcmp(newMode, 'ivp') ) % Going into IVP mode
    % IVP mode is almost identical to BVP, so do the same adjustments, and then
    % minor corrections at the end:
    handles = chebguiController.switchMode(handles, 'bvp');
    
    % Set the type of the problem
    handles.guifile.type = 'ivp';
    
    % Change heading for BC/IC input field:
    set(handles.panel_BCs, 'Title', 'Initial/final conditions')
    
    % Hide the discretization panel
    set(handles.panel_discretization, 'Visible', 'off');
    
    % Make the IVP solver panel visible
    set(handles.panel_IVPsolver, 'Visible', 'on')
    
    set(handles.button_bvp, 'Value', 0)
    set(handles.button_ivp, 'Value', 1)
    % Enable IVP solver option
    set(handles.menu_ivpSolver, 'Enable', 'on');
    
    % If we're solving in time-stepping mode, hide the initial guess box:
    if ( ~isempty(strfind(handles.guifile.options.ivpSolver, 'ode')) )
        set(handles.panel_initialGuess,'Visible','off');
    end
    
elseif ( strcmp(newMode, 'pde') ) % Going into PDE mode
    handles.guifile.type = 'pde';
    
    set(handles.button_bvp, 'Value', 0)
    set(handles.button_ivp, 'Value', 0)
    set(handles.button_pde, 'Value', 1)
    set(handles.button_eig, 'Value', 0)
    
    % Hide the iter_list box:
    handles = hideIterList(handles);
    
    % Show the initial guess panel, and ensure it has the correct title:
    set(handles.panel_initialGuess,'Visible','on','Title','Initial condition');
    
    set(handles.toggle_useLatest, 'Visible', 'off')
    set(handles.panel_DEs, 'Title', 'Differential equations')
    
    set(handles.input_GUESS, 'Enable', 'On')
    set(handles.toggle_useLatest, 'Value', 0)
    set(handles.toggle_useLatest, 'Enable', 'off')
    
    set(handles.panel_timedomain, 'Visible', 'on')
    
    set(handles.button_realplot, 'Visible', 'off')
    set(handles.button_imagplot, 'Visible', 'off')
    set(handles.popupmenu_bottomFig, 'Visible', 'off')

    set(handles.panel_eigopts, 'Visible', 'Off')
    
    % Switch from LBC and RBC fields to BC
    handles = toggleBCinput(handles, 'on', 'off');
    
    % Change the list of available options
    % Disable ODE menu options
    set(handles.menu_odedampednewton, 'Enable', 'Off')
    set(handles.menu_odeplotting, 'Enable', 'Off')    
    % Disable IVP solver option
    set(handles.menu_ivpSolver, 'Enable', 'off');    
    % Enable PDE menuoptions
    set(handles.menu_pdeSolver, 'Enable', 'On')
    set(handles.menu_pdeplotting, 'Enable', 'On')
    set(handles.menu_pdeholdplot, 'Enable', 'On')
    set(handles.menu_pdefix, 'Enable', 'On')
    set(handles.menu_fixN, 'Enable', 'On')
    
    % Clear the top figure
    chebguiController.initialiseFigureTop(handles)

    % Always clear the bottom figure
    chebguiController.initialiseFigureBottom(handles)
    
    % Hide the discretization and IVP solver panels
    set(handles.panel_discretization, 'Visible', 'off')
    set(handles.panel_IVPsolver, 'Visible', 'off')
    
    % Set correct header to top figure
    set(handles.panel_figSol, 'Title', 'Initial condition')
    
else % Going into EIG mode
    handles.guifile.type = 'eig';
    
    set(handles.button_bvp, 'Value', 0)
    set(handles.button_ivp, 'Value', 0)
    set(handles.button_pde, 'Value', 0)
    set(handles.button_eig, 'Value', 1)
    
    set(handles.toggle_useLatest, 'Visible', 'off')
    set(handles.panel_DEs, 'Title', 'Differential operator')
    set(handles.text_initial, 'String', 'Look for')
    set(handles.toggle_useLatest, 'Enable', 'off')
    
    set(handles.panel_timedomain, 'Visible', 'off')

    set(handles.button_realplot, 'Visible', 'on')
    set(handles.button_realplot, 'Value', 1)
    set(handles.button_imagplot, 'Visible', 'on')
    set(handles.button_imagplot, 'Value', 0)
    set(handles.popupmenu_bottomFig, 'Visible', 'off')

    set(handles.panel_eigopts,'Visible','On')
    if ( ~isempty(handles.guifile.sigma) )
        switch ( handles.guifile.sigma )
            case ''
                set(handles.popupmenu_sigma, 'Value', 1);
            case 'lm'
                set(handles.popupmenu_sigma, 'Value', 2);
            case 'sm'
                set(handles.popupmenu_sigma, 'Value', 3);
            case 'lr'
                set(handles.popupmenu_sigma, 'Value', 4);
            case 'sr'
                set(handles.popupmenu_sigma, 'Value', 5);
            case 'li'
                set(handles.popupmenu_sigma, 'Value', 6);
            case 'si'
                set(handles.popupmenu_sigma, 'Value', 7);
        end
    else
        set(handles.popupmenu_sigma, 'Value', 1)
    end

    if ( isfield(handles.guifile.options,'numeigs') ...
        && ~isempty(handles.guifile.options.numeigs) )
        set(handles.edit_eigN, 'String', handles.guifile.options.numeigs);
    else
        set(handles.edit_eigN, 'String', 6);
    end
    
    % Hide the initial guess panel, and ensure it has the correct title:
    set(handles.panel_initialGuess,'Visible','off');
    
    % Change heading for BC/IC input field:
    set(handles.panel_BCs, 'Title', 'Boundary conditions')

    % Hide the iter_list box:
    handles = hideIterList(handles);
    
    % Switch from LBC and RBC fields to BC
    handles = toggleBCinput(handles, 'off', 'on');
    
    % Change the list of available options.

    % Disable ODE options
    set(handles.menu_odedampednewton, 'Enable', 'Off')
    set(handles.menu_odeplotting, 'Enable', 'Off')
    
    % Disable IVP solver option
    set(handles.menu_ivpSolver, 'Enable', 'off');

    % Disable PDE options
    set(handles.menu_pdeSolver, 'Enable', 'Off')
    set(handles.menu_pdeplotting, 'Enable', 'Off')
    set(handles.menu_pdeholdplot, 'Enable', 'Off')
    set(handles.menu_pdefix, 'Enable', 'Off')
    set(handles.menu_fixN, 'Enable', 'Off')
    
    % Make the discretization panel visible
    set(handles.panel_discretization, 'Visible', 'on')
    
    % Hide the IVP solver panel
    set(handles.panel_IVPsolver, 'Visible', 'off')
    
    % Set correct header to top figure
    set(handles.panel_figSol, 'Title', 'Eigenvalues (imag vs real)')
    
    % Clear the top figure
    chebguiController.initialiseFigureTop(handles)
    
    % Always clear the bottom figure
    chebguiController.initialiseFigureBottom(handles)
    
    
end

end

function handles = toggleBCinput(handles, onOff1, onOff2)
set(handles.input_LBC, 'Visible', onOff1)
set(handles.panel_leftBCs, 'Visible', onOff1)
set(handles.input_RBC, 'Visible', onOff1)
set(handles.panel_rightBCs, 'Visible', onOff1)
set(handles.input_BC,  'Visible', onOff2)
set(handles.panel_BCs,  'Visible', onOff2)
end

function handles = hideIterList(handles)
set(handles.iter_list, 'Visible', 'off')
set(handles.iter_list, 'String', '')
set(handles.iter_list, 'Value', 0)
set(handles.iter_text, 'Visible', 'off')
set(handles.iter_text, 'String', '')
end


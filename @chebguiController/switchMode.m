function handles = switchMode(handles, newMode, callMode)
%SWITCHMODE   Go from one mode of CHEBGUI to another.
% Calling sequence:
%   HANDLES = SWITCHMODE(HANDLES, NEWMODE, CALLMODE)
% where
%   HANDLES:    A MATLAB handle object for the CHEBGUI figure.
%   NEWMODE:    Mode we want to switch to.
%   CALLMODE:   Indicates whether we are switching from a demo or not (in which
%               case, we migth want to show a plot of the initial guess of the
%               solution/initial condition).

% Copyright 2014 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% Only use callMode when calling from loadDemoMenu() - for demos with an initial
% guess / condition, we don't want to clear the figures.
if ( nargin == 2 )
    callMode = 'notDemo';
end

% Do a lot of disabling/enabling and hiding/showing objects on the CHEBGUI
% figure.
if ( strcmp(newMode, 'bvp') ) % Going into BVP mode
    handles.guifile.type = 'bvp';
    
    set(handles.button_ode, 'Value', 1)
    set(handles.button_pde, 'Value', 0)
    set(handles.button_eig, 'Value', 0)
    
    set(handles.text_DEs, 'String', 'Differential equation(s)')
    set(handles.input_GUESS, 'visible', 'On')
    set(handles.text_initial, 'String', 'Initial guess')
    set(handles.toggle_useLatest, 'Visible', 'on')
    
    set(handles.text_timedomain, 'Visible', 'off')
    set(handles.input_timedomain, 'Visible', 'off')
    
    set(handles.button_realplot, 'Visible', 'off')
    set(handles.button_imagplot, 'Visible', 'off')
    set(handles.popupmenu_bottomFig, 'Visible', 'on')
    set(handles.toggle_useLatest, 'Value', 0)
    
    if ( isempty(fieldnames(handles.latest)) )
        set(handles.toggle_useLatest, 'Enable', 'off')
    end
    
    set(handles.iter_list, 'Visible', 'off')
    set(handles.iter_list, 'String', '')
    set(handles.iter_list, 'Value', 0)
    set(handles.iter_text, 'Visible', 'off')
    set(handles.iter_text, 'String', '')
    
    set(handles.panel_eigopts, 'Visible', 'Off')
        
    % Switch from LBC and RBC fields to BC
    set(handles.input_LBC, 'Visible', 'Off')
    set(handles.text_LBCs, 'Visible', 'Off')
    set(handles.input_RBC, 'Visible', 'Off')
    set(handles.text_RBCs, 'Visible', 'Off')
    set(handles.input_BC, 'Visible', 'On')
    set(handles.text_BCs, 'Visible', 'On')
    
    % Change the list of available options.
    % Enable ODE menu options
    set(handles.menu_odedampednewton, 'Enable', 'On')
    set(handles.menu_odeplotting, 'Enable', 'On')
    % Disable PDE menu options
    set(handles.menu_pdeplotting, 'Enable', 'Off')
    set(handles.menu_pdeholdplot, 'Enable', 'Off')
    set(handles.menu_pdefix, 'Enable', 'Off')
    set(handles.menu_fixN, 'Enable', 'Off')
    
    % Clear the top figure
    if ( ~strcmp(callMode, 'demo') )
        chebguiController.initialiseFigureTop(handles)
    end
    % Always clear the bottom figure
    chebguiController.initialiseFigureBottom(handles)
    
    % Make the discretization panel visible
    set(handles.panel_discretization, 'Visible', 'on')
    
    % Note: The line below is a nice way to disable to objects, rather than hide
    % them.
    % Enable all childrens of discretisation panel
    %     set(findall(handles.panel_discretization,  '-property',  'enable'),  ...
    %         'visible',  'on')
    
elseif ( strcmp(newMode, 'pde') ) % Going into PDE mode
    handles.guifile.type = 'pde';
    
    set(handles.button_ode, 'Value', 0)
    set(handles.button_pde, 'Value', 1)
    set(handles.button_eig, 'Value', 0)
    set(handles.toggle_useLatest, 'Visible', 'off')
    set(handles.iter_list, 'Visible', 'off')
    set(handles.iter_text, 'Visible', 'off')
    set(handles.text_initial, 'String', 'Initial condition')
    set(handles.text_DEs, 'String', 'Differential equation(s)')
    
    set(handles.input_GUESS, 'Visible', 'On')
    set(handles.input_GUESS, 'Enable', 'On')
    set(handles.toggle_useLatest, 'Value', 0)
    set(handles.toggle_useLatest, 'Enable', 'off')
    
    set(handles.text_timedomain, 'Visible', 'on')
    set(handles.input_timedomain, 'Visible', 'on')
    
    set(handles.button_realplot, 'Visible', 'off')
    set(handles.button_imagplot, 'Visible', 'off')
    set(handles.popupmenu_bottomFig, 'Visible', 'off')

    set(handles.panel_eigopts, 'Visible', 'Off')
    
    % Switch from BC to LBC and RBC fields
    set(handles.input_LBC, 'Visible', 'On')
    set(handles.text_LBCs, 'Visible', 'On')
    set(handles.input_RBC, 'Visible', 'On')
    set(handles.text_RBCs, 'Visible', 'On')
    set(handles.input_BC, 'Visible', 'Off')
    set(handles.text_BCs, 'Visible', 'Off')
    
    % Change the list of available options
    % Disable ODE menu options
    set(handles.menu_odedampednewton, 'Enable', 'Off')
    set(handles.menu_odeplotting, 'Enable', 'Off')
    % Enable PDE menuoptions
    set(handles.menu_pdeplotting, 'Enable', 'On')
    set(handles.menu_pdeholdplot, 'Enable', 'On')
    set(handles.menu_pdefix, 'Enable', 'On')
    set(handles.menu_fixN, 'Enable', 'On')
    
    % Clear the top figure
    if ( ~strcmp(callMode, 'demo') )
        chebguiController.initialiseFigureTop(handles)
    end
    % Always clear the bottom figure
    chebguiController.initialiseFigureBottom(handles)
    
    % Hide the discretization panel
    set(handles.panel_discretization, 'Visible', 'off')
    
    
else % Going into EIG mode
    handles.guifile.type = 'eig';
    
    set(handles.button_ode, 'Value', 0)
    set(handles.button_pde, 'Value', 0)
    set(handles.button_eig, 'Value', 1)
    
    set(handles.toggle_useLatest, 'Visible', 'off')
    set(handles.text_DEs, 'String', 'Differential operator')
    % set(handles.iter_list, 'Visible', 'on')
    % set(handles.iter_text, 'Visible', 'on')
    set(handles.text_initial, 'String', 'Look for')
    set(handles.input_GUESS, 'visible', 'Off')
    set(handles.toggle_useLatest, 'Enable', 'off')
    
    set(handles.text_timedomain, 'Visible', 'off')
    set(handles.input_timedomain, 'Visible', 'off')

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
    
    set(handles.iter_list, 'Visible', 'off')
    set(handles.iter_list, 'String', '')
    set(handles.iter_list, 'Value', 0)
    set(handles.iter_text, 'Visible', 'off')
    set(handles.iter_text, 'String', '')
    
    % Switch from LBC and RBC fields to BC
    set(handles.input_LBC, 'Visible', 'Off')
    set(handles.text_LBCs, 'Visible', 'Off')
    set(handles.input_RBC, 'Visible', 'Off')
    set(handles.text_RBCs, 'Visible', 'Off')
    set(handles.input_BC, 'Visible', 'On')
    set(handles.text_BCs, 'Visible', 'On')
    
    % Change the list of available options.

    % Disable ODE options
    set(handles.menu_odedampednewton, 'Enable', 'Off')
    set(handles.menu_odeplotting, 'Enable', 'Off')

    % Disable PDE options
    set(handles.menu_pdeplotting, 'Enable', 'Off')
    set(handles.menu_pdeholdplot, 'Enable', 'Off')
    set(handles.menu_pdefix, 'Enable', 'Off')
    set(handles.menu_fixN, 'Enable', 'Off')
    
    
    % Make the discretization panel visible
    set(handles.panel_discretization, 'Visible', 'on')
    
    % Clear the top figure
    if ( ~strcmp(callMode, 'demo') )
        chebguiController.initialiseFigureTop(handles)
    end
    % Always clear the bottom figure
    chebguiController.initialiseFigureBottom(handles)
    
    
end

end

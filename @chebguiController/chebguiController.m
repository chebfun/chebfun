classdef chebguiController
%CHEBGUICONTROLLER   Control the layout of CHEBGUI.
%   This class is not intended to be called directly by the end user.
%
% See also CHEBGUI.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEVELOPER NOTE
%   The CHEBGUICONTROLLER class implements a number of methods, used to control
%   the look of CHEBGUI. In v4, this functionality used to live in the @chebgui
%   folder, but to increase modularity, it has been spun off to its own class.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CLASS METHODS:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods ( Access = public, Static = false )
        
        function A = chebguiController(varargin)
            % We never construct objects of this type, so the constructor is
            % trivial
        end
        
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% STATIC METHODS:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods ( Access = public, Static = true )
        
        % Deal with what happens when BCs get changed.
        handles = callbackBCs(handles, inputString, type)
        
        % Change font size in CHEBGUI
        handles = changeFontsize(handles, change)
        
        % Clear everything in the CHEBGUI window
        handles = clear(handles)
        
        % Draw the Chebfun logo on the GUI:
        handles = drawLogo(handles)

        function initialiseFigureBottom(handles)
            %INITIALISEFIGUREBOTTOM    Reset bottom figure of CHEBGUI.

            cla(handles.fig_norm, 'reset');
            axes(handles.fig_norm)
            box on
            set(handles.fig_norm, 'fontsize', handles.fontsizePanels);
            set(handles.panel_figNorm, 'title', ' ');
        end
        
        function initialiseFigureTop(handles)
            %INITIALISEFIGURETOP    Reset top figure of CHEBGUI.
            
            cla(handles.fig_sol, 'reset');
            axes(handles.fig_sol)
            box on
            set(handles.fig_sol, 'fontsize', handles.fontsizePanels);
            set(handles.panel_figNorm, 'title', 'Solution');

        end

        function initialiseFigures(handles)
            %INITIALISEFIGURES    Reset top and bottom figures of CHEBGUI.
            
            chebguiController.initialiseFigureBottom(handles);
            chebguiController.initialiseFigureTop(handles);
        end
        
        function handles = initialiseMenus(handles)
            %INITALISEMENUS    Create menu items for CHEBGUI
            
            % Don't want to create the menus if they have already been created
            % (e.g. if we call CHEBGUI again from the command line when it's
            % already open):
            if ( isfield(handles,'demosLoaded') )
                return
            end
            
            % For ODE problems
            handles.menu_bvps = uimenu(handles.menu_demos, ...
                'label', 'ODE - Scalar BVPs');
            handles.menu_systems = uimenu(handles.menu_demos, ...
                'label', 'ODE - Coupled BVPs');
            handles.menu_ivps = uimenu(handles.menu_demos, ...
                'label', 'ODE - Scalar IVPs', 'separator', 'on');
            handles.menu_IVPsystems = uimenu(handles.menu_demos, ...
                'label', 'ODE - Coupled IVPs');

            % For EIG problems
            handles.menu_eigsscalar = uimenu(handles.menu_demos, ...
                'label', 'EIG - Scalar','separator','on');
            handles.menu_eigssystem = uimenu(handles.menu_demos, ...
                'label', 'EIG - Systems');
            
            % For PDE problems
            handles.menu_pdesingle = uimenu(handles.menu_demos, ...
                'label', 'PDE - Scalar','separator','on');
            handles.menu_pdesystems = uimenu(handles.menu_demos, ...
                'label', 'PDE - Systems');
                        
        end
                
        % Initialize fonts of the CHEBGUI window
        handles = initalizeFields(handles)
        
        % Plot eigenmodes in the GUI
        plotEigenmodes(handles, selection, h1, h2)
        
        % Populate the fields of the CHEBGUI figure.
        populate(hObject, handles, guifile)
        
        % Load the menu in CHEBGUI of demos.
        loadDemoMenu(handles)
        
        % Remove tabs from input
        str = removeTabs(str)
        
        % Set-up the panels on CHEBGUI
        handles = setupPanels(handles);

        % Set-up the figures panel on CHEBGUI
        handles = setupPanelFigures(handles);
        
        % Set-up the input panel on CHEBGUI
        handles = setupPanelInput(handles);

        % Set-up the discretization option panel on CHEBGUI
        handles = setupPanelDiscretization(handles);
        
        % Set-up the IVP solver option panel on CHEBGUI
        handles = setupPanelIVPsolver(handles);
        
        % Set-up the panels on CHEBGUI
        handles = setupPanelType(handles);
        
        % Switch between different modes in CHEBGUI (BVP, EIG or PDE).
        handles = switchMode(handles, newMode)
        
    end
    
end
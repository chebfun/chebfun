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
        
        % Clear everything in the CHEBGUI window
        handles = clear(handles)

        function initialiseFigureBottom(handles)
            %INITIALISEFIGUREBOTTOM    Reset bottom figure of CHEBGUI.

            cla(handles.fig_norm, 'reset');
            axes(handles.fig_norm)
            box on
        end
        
        function initialiseFigureTop(handles)
            %INITIALISEFIGURETOP    Reset top figure of CHEBGUI.

            cla(handles.fig_sol, 'reset');
            axes(handles.fig_sol)
            title('Solution')
            box on

        end

        function initialiseFigures(handles)
            %INITIALISEFIGURES    Reset top and bottom figures of CHEBGUI.
            
            chebguiController.initialiseFigureBottom(handles);
            chebguiController.initialiseFigureTop(handles);
        end

        
        
        % Initialize fonts of the CHEBGUI window
        handles = initalizeFields(handles)
        
        % Plot eigenmodes in the GUI
        plotEigenmodes(handles, selection, h1, h2)
        
        % Populate the fields of the CHEBGUI figure.
        initSuccess = populate(handles, guifile)
        
        % Load the menu in CHEBGUI of demos.
        loadDemoMenu(handles)
        
        % Switch between different modes in CHEBGUI (BVP, EIG or PDE).
        handles = switchMode(handles, newMode, callMode)
        
    end
    
end
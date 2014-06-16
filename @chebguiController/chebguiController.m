classdef chebguiController
%CHEBGUICONTROLLER   Control the layout of CHEBGUI.
%   This class is not intended to be called directly by the end user.
%
%   See also CHEBGUI.

% Developers note:
%   The CHEBGUICONTROLLER class implements a number of methods, used to control
%   the look of CHEBGUI. In v4, this functionality used to live in the @chebgui
%   folder, but to increase modularity, it has been spun off to its own class.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.
    
    %% Non-Static methods:
    methods ( Access = public, Static = false )
        
        function A = chebguiController(varargin)
            % We never construct objects of this type, so the constructor is
            % trivial
        end
        
    end
    
    %% Static methods:
    methods ( Access = public, Static = true )
        
        % Deal with what happens when BCs get changed.
        handles = callbackBCs(handles, inputString, type)
        
        % Clear everything in the CHEBGUI window
        handles = clear(handles)
        
        % Load the menu in CHEBGUI of demos.
        loadDemoMenu(handles)
        
        % Switch between different modes in CHEBGUI (BVP, EIG or PDE).
        handles = switchMode(handles, newMode, callMode)
        
        % Plot eigenmodes in the GUI
        plotEigenmodes(handles, selection, h1, h2)
       
        % Populate the fields of the CHEBGUI figure.
        initSuccess = populate(handles, guifile)
        
    end
    
end
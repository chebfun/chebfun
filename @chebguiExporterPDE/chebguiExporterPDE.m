classdef chebguiExporterPDE < chebguiExporter
%CHEBGUIEXPORTERPDE   Export a PDE from CHEBGUI.
%   This is a an concrete implementation of the class CHEBGUIEXPORTER, which
%   exports PDEs from CHEBGUI to .m-files, to the workspace, or to a .chebgui
%   file. It is not intended to be called directly by the end user.
%
%   See also CHEBGUI, CHEBGUIEXPORTER.

% Developers note:
%   The CHEBGUICONTROLLER class defines a number of abstract methods, used to
%   export problems from CHEBGUI. In v4, this functionality used to live in the
%   @chebgui folder, but to increase modularity, it has been spun off to its own
%   class.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.
        
    methods (Access = public)
        
        function A = chebguiExporterPDE(varargin)
            % Set default file name:
            A.defaultFileName = 'pde.m';
        end
        
    end
    
    
    methods ( Static = true )
        
        function e = make(varargin)
            % Factory method.
            e = chebguiExporterPDE(varargin{:}); 
        end
        
    end
end
classdef chebguiExporterPDE < chebguiExporter
    %CHEBGUIEXPORTERPDE   Export a PDE from CHEBGUI.
    %   This is a an concrete implementation of the class CHEBGUIEXPORTER, which
    %   exports PDEs from CHEBGUI to .m-files, to the workspace, or to a
    %   .chebgui file. It is not intended to be called directly by the end user.
    %
    %   See also CHEBGUI, CHEBGUIEXPORTER.
    
    % Copyright 2014 by The University of Oxford and The Chebfun Developers.
    % See http://www.chebfun.org/ for Chebfun information.
    
    properties
    
        % The default file name when exporting to an .m-file:
        defaultFileName = 'pde.m';
        
        % Description for printing to .m files:
        description = 'a partial differential equation.';
    
    end
    
    methods (Access = public)
        
        function A = chebguiExporterPDE(varargin)
            % Do nothing!
        end
        
    end
    
    
    methods ( Static = true )
        
        % Extract information from the CHEBGUI object to a struct
        expInfo = exportInfo(guifile)
        
        % Print problem description:
        printDescription(fid, expInfo)
        
        % Print options for solving the problem:
        printOptions(fid, expInfo)
        
        % Print lines for setting up the problem:
        printSetup(fid, expInfo, guifile)
        
        % Print the actual lines for calling the solver method:
        printSolver(fid, expInfo)
        
        % Print steps taken after the solver finishes:
        printPostSolver(fid, expInfo)
        
    end
end
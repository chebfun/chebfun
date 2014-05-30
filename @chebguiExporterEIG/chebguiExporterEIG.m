classdef chebguiExporterEIG < chebguiExporter
    %CHEBGUIEXPORTEREIG   Export an EIG problem from CHEBGUI.
    %   This is a an concrete implementation of the class CHEBGUIEXPORTER, which
    %   exports EIG problems from CHEBGUI to .m-files, to the workspace, or to a
    %   .chebgui file. It is not intended to be called directly by the end user.
    %
    %   See also CHEBGUI, CHEBGUIEXPORTER.
    
    % Copyright 2014 by The University of Oxford and The Chebfun Developers.
    % See http://www.chebfun.org/ for Chebfun information.
    
    properties
        
        % The default file name when exporting to an .m-file:
        defaultFileName = 'bvpeig.m';
        
        % Description for printing to .m files:
        description = 'an eigenvalue problem';
        
    end
    
    methods (Access = public)
        
        function A = chebguiExporterEIG(varargin)
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
        
        function toWorkspaceSolutionOnly(handles)
            varnames = handles.varnames;
            lambdaName = handles.eigVarName;
            
            if ( iscell(lambdaName) )
                lambdaName = lambdaName{:};
            end
            
            nv = numel(varnames);
            d = handles.latest.solution;
            V = handles.latest.solutionT;
            if ( ~iscell(V) )
                V = {V};
            end
            
            for k = 1:nv
                assignin('base', varnames{k}, V{k});
            end
            assignin('base', lambdaName, d);
            evalin('base', lambdaName);
        end
        
        function toMat(handles)
            D = diag(handles.latest.solution); %#ok<NASGU>
            V = handles.latest.solutionT;  %#ok<NASGU>
            uisave({'D', 'V'}, 'bvpeig');
        end
    end
end
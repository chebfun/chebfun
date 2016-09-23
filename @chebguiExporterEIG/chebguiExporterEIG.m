classdef chebguiExporterEIG < chebguiExporter
%CHEBGUIEXPORTEREIG   Export an EIG problem from CHEBGUI.
%   This is a an concrete implementation of the class CHEBGUIEXPORTER, which
%   exports EIG problems from CHEBGUI to .m-files, to the workspace, or to a
%   .chebgui file. It is not intended to be called directly by the end user.
%
% See also CHEBGUI, CHEBGUIEXPORTER.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CLASS PROPERTIES:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    properties ( Access = public )
        
        % The default file name when exporting to an .m-file:
        defaultFileName = 'chebeig.m';
        
        % Description for printing to .m files:
        description = 'an eigenvalue problem';
        
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CLASS CONSTRUCTOR:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods ( Access = public, Static = false )
        
        function A = chebguiExporterEIG(varargin)
            % Do nothing!
        end
        
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% STATIC METHODS:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods ( Access = public, Static = true )
        
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
            %TOWORKSPACESOLUTIONONLY    Export the solution to the workspace
            
            % Get the variable names involved in the problem:
            varnames = handles.varnames;
            lambdaName = handles.eigVarName;
            
            if ( iscell(lambdaName) )
                lambdaName = lambdaName{:};
            end
                        
            % Extract the latest solution
            d = handles.latest.solution;
            U = handles.latest.solutionT;
            
            % Export to workspace
            assignin('base', lambdaName, d);
            assignin('base', 'U', U);
            % Print in the command window the name of the variables
            evalin('base', lambdaName);
            evalin('base', 'U');
        end
        
        function toMat(handles)
            %TOMAT  Export from the CHEBGUI figure to a .mat-file.
            
            % Obtain the variables involved in the problem.
            D = diag(handles.latest.solution); %#ok<NASGU>
            V = handles.latest.solutionT;  %#ok<NASGU>
            
            % Ask user where to save:
            uisave({'D', 'V'}, 'bvpeig');
        end
        
        function toWorkspace(handles)
            %TOWORKSPACE    Export from the CHEBGUI figure to the workspace.
            
            % Setup dialog for asking the user for variable names to be used.
            prompt = {'Eigenvalues', 'Eigenmodes'};
            name = 'Export to workspace';
            defaultAnswer ={'lambda', 'U'};
            numlines = 1;
            options.Resize ='on';
            options.WindowStyle ='modal';
            
            % Ask user what variable names to use:
            answer = inputdlg(prompt, name, numlines, defaultAnswer, options);
            
            % Assign to workspace:
            if ( ~isempty(answer) )
                assignin('base', answer{1}, handles.latest.solution);
                assignin('base', answer{2}, handles.latest.solutionT);
                % Print in the command window the name of the variables
                evalin('base', answer{1})
                evalin('base', answer{2});
            end
        end
        
    end
    
end

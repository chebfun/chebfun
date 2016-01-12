classdef chebguiExporterPDE < chebguiExporter
%CHEBGUIEXPORTERPDE   Export a PDE from CHEBGUI.
%   This is a an concrete implementation of the class CHEBGUIEXPORTER, which
%   exports PDEs from CHEBGUI to .m-files, to the workspace, or to a .chebgui
%   file. It is not intended to be called directly by the end user.
%
% See also CHEBGUI, CHEBGUIEXPORTER.

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CLASS PROPERTIES:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    properties ( Access = public )
        
        % The default file name when exporting to an .m-file:
        defaultFileName = 'chebpde.m';
        
        % Description for printing to .m files:
        description = 'a partial differential equation';
        
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CLASS CONSTRUCTOR:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods ( Access = public, Static = false )
        
        function A = chebguiExporterPDE(varargin)
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
            %TOWORKSPACESOLUTIONONLY   Export the solution to the workspace.
            %
            % Note: this method exports fewer objects to the workspace than the
            % toWorkspace() method.
            
            % Obtain the names of the variables involved:
            varnames = handles.varnames;
            nv = numel(varnames);
            
            % Obtain the latest solution:
            sol = handles.latest.solution;
            if ( ~iscell(sol) )
                sol = {sol};
            end
            
            % Assign to workspace:
            for k = 1:nv
                assignin('base', varnames{k}, sol{k});
                evalin('base', varnames{k});
            end
        end
        
        function toMat(handles)
            %TOMAT   Export from the CHEBGUI figure to a .mat-file.
            
            % Obtain the variable names involved in the problem.
            varnames = handles.varnames;
            
            % Obtain the latest solution:
            sol = handles.latest.solution;
            if ( ~iscell(sol) )
                sol = {sol};
            end
            
            % Assign variables before exporting:
            for k = 1:numel(varnames);
                eval([varnames{k} ' = sol{k};']);
            end
            t = handles.latest.solutionT;  %#ok<NASGU>
            
            % Ask user where to save the .mat file:
            uisave([varnames', 't'], 'pde');
        end
        
        function toWorkspace(handles)
            %TOWORKSPACE   Export from the CHEBGUI figure to the workspace.
            
            % Setup dialog for asking the user for variable names to be used.
            varnames = handles.varnames;
            nv = numel(varnames);
            if ( nv == 1 )
                prompt = {'Solution', 'Time domain', 'Final solution'};
                defaultAnswer=  [varnames,'t', [varnames{1} '_final']];
            else
                sol = {};
                solfinal = {};
                finalsol = {};
                for k = 1:nv
                   BVP % Show the expected names in the dialog:
                    sol = [sol ['Solution ' varnames{k}]];
                    finalsol = [finalsol ['Final ' varnames{k}]];
                    solfinal = [solfinal [varnames{k} '_final']];
                end
                prompt = [sol, 'Time domain', finalsol];
                defaultAnswer = [varnames', 't', solfinal];
            end
            
            name = 'Export to workspace';
            numlines = 1;
            options.Resize ='on';
            options.WindowStyle ='modal';
            
            % Ask user what variable names he/she wants to use:
            answer = inputdlg(prompt, name, numlines, defaultAnswer, options);
            
            % Obtain the latest solution:
            sol = handles.latest.solution;
            if ( ~iscell(sol) )
                sol = {sol};
            end
            
            % Assign to the workspace!
            if ( ~isempty(answer) )
                for k = 1:nv
                    assignin('base', answer{k}, sol{k});
                end
                assignin('base', answer{nv+1}, handles.latest.solutionT);
                for k = 1:nv
                    assignin('base', answer{nv+1+k}, sol{k}(:,end));
                end
            end
            
        end
        
    end
    
end

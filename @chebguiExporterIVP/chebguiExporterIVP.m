classdef chebguiExporterIVP < chebguiExporter
%CHEBGUIEXPORTERIVP   Export an IVP from CHEBGUI.
%   This is a an concrete implementation of the class CHEBGUIEXPORTER, which
%   exports IVPs from CHEBGUI to .m-files, to the workspace, or to a .chebgui
%   file. It is not intended to be called directly by the end user.
%
% See also CHEBGUI, CHEBGUIEXPORTER.

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CLASS PROPERTIES:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    properties
        
        % The default file name when exporting to an .m-file:
        defaultFileName = 'chebivp.m';
        
        % Description for printing to .m files:
        description = 'an initial-value problem';
        
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CLASS CONSTRUCTOR:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods ( Access = public, Static = false )
        
        function A = chebguiExporterIVP(varargin)
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
        %TOWORKSPACESOLUTIONONLY   Export the solution to the workspace
        %
        % Note: this method exports fewer objects to the workspace than the
        % toWorkspace() method.
            
            % Obtain the variable names:
            varnames = handles.varnames;
            nv = numel(varnames);
            
            % Obtain the latest solution:
            sol = handles.latest.solution;
            
            % Export to the workspace!
            for k = 1:nv
                assignin('base', varnames{k}, sol(:,k));
                evalin('base', varnames{k});
            end
        end
        
        function toMat(handles)
        %TOMAT   Export from the CHEBGUI figure to a .mat-file.
            
            % Obtain the variable names involved in the problem.
            varnames = handles.varnames;
            for k = 1:numel(varnames);
                eval([varnames{k} ' = handles.latest.solution(:,k);']);
            end
            
            % Extract the norm of the updates, the CHEBOP and all options.
            N = handles.latest.chebop;  %#ok<NASGU>
            options = handles.latest.options;  %#ok<NASGU>
            
            % Show the user a figure for selecting where to save.
            uisave([varnames', 'N', 'options'], 'ivp');
        end
        
        function toWorkspace(handles)
        %TOWORKSPACE   Export from the CHEBGUI figure to the workspace.
            
            % Setup dialog for asking the user for variable names to be used.
            numlines = 1;
            options.Resize ='on';
            options.WindowStyle ='modal';
            
            % Need slightly different dialogs depending on whether we were
            % solving scalar problem or a coupled system.
            varnames = handles.varnames;
            nv = numel(varnames);
            if ( nv == 1 )
                prompt = {'Differential operator', 'Solution:', 'Options'};
            else
                prompt = ['Differential operator', varnames.', 'Options'];
            end
            
            name = 'Export to workspace';
            
            defaultAnswer = ['N', varnames', 'options'];
            
            % Show the user a dialog.
            answer = inputdlg(prompt, name, numlines, defaultAnswer, options);
            
            % Obtain the latest solution:
            sol = handles.latest.solution;
            
            % Export to the workspace!
            if ( ~isempty(answer) )
                assignin('base', answer{1}, handles.latest.chebop);
                for k = 1:nv
                    assignin('base', answer{k+1}, sol(:,k));
                end
                assignin('base', answer{nv+2}, handles.latest.options);
            end
        end
        
    end
    
end

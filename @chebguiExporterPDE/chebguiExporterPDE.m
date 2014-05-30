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
        
        function toWorkspaceSolutionOnly(handles)
            varnames = handles.varnames;
            nv = numel(varnames);
            
            sol = handles.latest.solution;
            if ( ~iscell(sol) )
                sol = {sol};
            end
            
            for k = 1:nv
                assignin('base', varnames{k}, sol{k});
                evalin('base', varnames{k});
            end
        end
        
        function toMat(handles)
            varnames = handles.varnames;
            
            sol = handles.latest.solution;
            if ( ~iscell(sol) )
                sol = {sol};
            end
            
            for k = 1:numel(varnames);
                eval([varnames{k} ' = sol{k};']);
            end
            t = handles.latest.solutionT;  %#ok<NASGU>
            uisave([varnames', 't'], 'pde');
        end
        
        function toWorkspace(handles)
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
            
            answer = inputdlg(prompt, name, numlines, defaultAnswer, options);
            
            sol = handles.latest.solution;
            if ( ~iscell(sol) )
                sol = {sol};
            end
            
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
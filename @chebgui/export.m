function export(guifile, handles, exportType)

% TODO:  Documentation.

% Copyright 2014 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

problemType = guifile.type;

% handles and exportType are empty if user calls export from command window.
% Always export to .m file in that case.
if ( nargin == 1 )
    exportType = '.m';
end

% Exporting a BVP
if ( strcmp(problemType, 'bvp') )
    switch exportType
        case 'Workspace'
            numlines = 1;
            options.Resize ='on';
            options.WindowStyle ='modal';
                
            varnames = handles.varnames;
            nv = numel(varnames);
            if ( nv == 1 )
                prompt = {'Differential operator', 'Solution:',...
                'Vector with norm of updates', 'Options'};
            else
                prompt = ['Differential operator', varnames.',...
                'Vector with norm of updates', 'Options'];
            end
                    
            name = 'Export to workspace';
            
            defaultAnswer = ['N', varnames', 'normVec', 'options'];

            answer = inputdlg(prompt, name, numlines, defaultAnswer, options);
            
            sol = handles.latest.solution;
            if ( ~isempty(answer) )
                assignin('base', answer{1}, handles.latest.chebop);
                for k = 1:nv
                    assignin('base', answer{k+1}, sol(:,k));
                end
                assignin('base', answer{nv+2}, handles.latest.norms);
                assignin('base', answer{nv+3}, handles.latest.options);
            end
            
        case 'Cancel'
            return;
    end
    
% Exporting a PDE
elseif ( strcmp(problemType,'pde') )
    switch exportType
        case 'Workspace'           
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
            
        case 'Cancel'
            return;
    end

% Exporting an EIG problem.
else
    switch exportType
        case 'Workspace'           
            prompt = {'Eigenvalues', 'Eigenmodes'};
            name = 'Export to workspace';
            defaultAnswer ={'D', 'V'};
            numlines = 1;
            options.Resize ='on';
            options.WindowStyle ='modal';
            
            answer = inputdlg(prompt, name, numlines, defaultAnswer, options);
            
            if ( ~isempty(answer) )
                assignin('base', answer{1}, diag(handles.latest.solution));
                assignin('base', answer{2}, handles.latest.solutionT);
            end

        case 'Cancel'
            return;
    end
end

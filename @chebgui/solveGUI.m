function handles = solveGUI(guifile, handles)
%SOLVEGUI    Called when a user hits the solve button of the Chebfun GUI.
%
% Calling sequence:
%
%   HANDLES = SOLVEGUI(GUIFILE, HANDLES)
%
% where
%   
%   HANDLES:    A MATLAB handle object to the CHEBGUIWINDOW figure.
%   GUIFILE:    A CHEBGUI object.
%
% See also: chebgui/solveGUIbvp, chebgui/solveGUIeig, chebgui/solveGUIivp,
%           chebgui/solveGUIpde.

% Copyright 2014 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% Check whether some input is missing
if ( isempty(guifile.domain) )
    error('CHEBFUN:CHEBGUI:solveGUI:emptyDomain', ...
        'The domain must be defined.');
end

if ( isempty(guifile.DE) )
    error('CHEBFUN:CHEBGUI:solveGUI:emptyDE', ...
        'The differential equation can not be empty.');
end

if ( isempty(guifile.LBC) && isempty(guifile.RBC) && isempty(guifile.BC) )
    error('CHEBFUN:CHEBGUI:solveGUI:emptyBC', ...
        'Boundary conditions must be defined.');
end

if ( strcmp(get(handles.button_solve, 'string'), 'Solve') )   % In solve mode
    
    % Some basic checking of the inputs.
    dom = guifile.domain;
    a = dom(1);
    b = dom(end);
    if ( b <= a )
        error('CHEBFUN:CHEBGUI:solveGUI:incorrectDomain', ...
            'Error in constructing domain. %s is not valid.',dom);
    end

    % Check that we have tolerance information:
    tol = guifile.tol;
    if ( ~isempty(tol) )
        tolnum = str2double(tol);
        if ( isnan(tolnum) || isinf(tolnum) || isempty(tolnum) )
            error('CHEBFUN:CHEBGUI:solveGUI:invalidTolerance', ...
                'Invalid tolerance, ''%s''.', tol);
        end
    end   

    if ( strcmpi(guifile.type,'pde') )
        % If we're in PDE mode, need to have time information.
        tt = str2num(guifile.timedomain);
        if ( isempty(tt) )
            error('CHEBFUN:CHEBGUI:solveGUI:incorrectTimeInterval', ...
                'Error in constructing time interval.');
        end
        if ( isempty(guifile.init) )
            error('CHEBFUN:CHEBGUI:solveGUI:emptyInitialCondition', ...
                'Initial condition is empty.');
        end
        if ( str2double(tol) < 1e-6 )
            tolchk = questdlg(['WARNING: PDE solves in chebgui are limited ' ...
                'to a tolerance of 1e-6'], ...
                'WARNING','Continue', 'Cancel','Continue');
            if ( strcmp(tolchk,'Continue') )
                tol = '1e-6';
                guifile.tol = tol;
                handles.guifile.tol = tol;
            else
                return
            end
        end
    end
    
    % Disable buttons, figures, etc.
    set(handles.toggle_useLatest, 'Enable','off');
    set(handles.button_exportsoln, 'Enable','off');
    set(handles.button_figsol, 'Enable','off');

    % STOP and PAUSE don't work in EIGS mode.
    if ( ~get(handles.button_eig, 'Value') )
        % Pause button
        set(handles.button_clear, 'String','Pause');
        set(handles.button_clear, 'BackgroundColor',[255 179 0]/256);
        % Stop button
        set(handles.button_solve, 'String','Stop');
        set(handles.button_solve, 'BackgroundColor',[214 80 80]/256);
    end
    
    % Update the figure:
    drawnow
    set(handles.menu_demos, 'Enable','off');

    % What discretization do we want to use?
    if ( get(handles.button_collocation, 'Value') )
        guifile.options.discretization = 'collocation';
    else
        guifile.options.discretization = 'ultraspherical';
    end
    
    % Call the private method solveGUIbvp, ivp, pde, or eig, which do the work.
    try
        if ( strcmpi(handles.guifile.type, 'bvp') )
            handles = solveGUIbvp(guifile, handles);
        elseif ( strcmpi(handles.guifile.type, 'ivp') )
            handles = solveGUIivp(guifile, handles);
        elseif ( strcmpi(handles.guifile.type, 'pde') )
            handles = solveGUIpde(guifile, handles);
        else
            handles = solveGUIeig(guifile, handles);            
        end
        handles.hasSolution = 1;
    catch ME
        Mstruct = struct('identifier', ME.identifier, 'message', ME.message, ...
            'stack', ME.stack);
        MEID = ME.identifier;

        % For specific common errors, update the error information before
        % rethrowing an error.
        if ( strcmp(MEID,'LINOP:mldivide:NoConverge') )
            Mstruct.identifier = ['Chebgui:' Mstruct.identifier];
            Mstruct.message = [Mstruct.message, ...
                'See "help cheboppref" for details on how to increase ' ...
                'number of points.'];
        elseif ( ~isempty(strfind(MEID,'Parse:')) ...
                || ~isempty(strfind(MEID,'LINOP:')) ...
                ||~isempty(strfind(MEID,'Lexer:')) ...
                || ~isempty(strfind(MEID,'Chebgui:')) )
            Mstruct.identifier = ['Chebgui:' Mstruct.identifier];
        elseif ( strcmp(MEID,'CHEBOP:solve:findguess:DivisionByZeroChebfun') )
            Mstruct.identifier = ['Chebgui:' Mstruct.identifier];
            Mstruct.message =  ['Error in constructing initial guess. The ' ...
                'zero function on the domain is not a permitted initial ' ...
                'guess as it causes division by zero. Please assign an ' ...
                'initial guess using the initial guess field.'];
        else
            handles.hasSolution = 0;
            rethrow(ME);
        end
        error(Mstruct)
    end

    % Once we have finished solving, update parts of the GUI to reflect it.
    resetComponents(handles);
    set(handles.toggle_useLatest, 'Enable', 'on');
    set(handles.button_exportsoln, 'Enable', 'on');
else   % In stop mode
    set(handles.button_clear, 'String', 'Clear all');
    set(handles.button_clear, 'BackgroundColor', ...
        get(handles.button_export,'BackgroundColor'));
    set(handles.button_solve, 'String','Solve');
    set(handles.button_solve, 'BackgroundColor', [43 129 86]/256);
    set(handles.menu_demos, 'Enable','on');
    set(handles.button_exportsoln, 'Enable', 'off');
    handles.hasSolution = 1;
    drawnow
end

end

function resetComponents(handles)
%RESETCOMPONENTS    Update GUI to reflect that we have solved a problem

% Enable buttons, figures, etc. Set button to 'solve' again
set(handles.button_solve, 'String', 'Solve');
set(handles.button_solve, 'BackgroundColor', [43 129 86]/256);
set(handles.button_clear, 'String', 'Clear all');
set(handles.button_clear, 'BackgroundColor', ...
    get(handles.button_export, 'BackgroundColor'));
set(handles.button_figsol, 'Enable', 'on');
set(handles.button_exportsoln, 'Enable', 'off');
set(handles.menu_demos, 'Enable', 'on');

end

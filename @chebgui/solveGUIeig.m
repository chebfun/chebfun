function varargout = solveGUIeig(guifile, handles)
%SOLVEGUIEIG   Solve a eigenvalue problem, specified by a CHEBGUI object.
%
% Calling sequence:
%
%   VARARGOUT = SOLVEGUIEIG(GUIFILE, HANDLES)
%
% where
%   
%   GUIFILE:    A CHEBGUI object, describing the problem.
%   HANDLES:    A MATLAB handle to the chebguiwindow figure.
%
% If the method is called by pressing the 'Solve' button on the GUI,
%   VARARGOUT{1}:   Will be a MATLAB handle to the chebguiwindow figure, which
%                   has been updated to contain the solution and other useful
%                   results for the problem solved.
%
% If the method is called by calling the command explicitly with a CHEBGUI
% object (e.g. [V, D] = SOLVEGUIEIG(GUIFILE) from the command line),
%   VARARGOUT{1}:   A diagonal matrix containing the eigenvalues.
%   VARARGOUT{2}:   A CHEBMATRIX of the eigenfunctions.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Handles will be an empty variable if we are solving without using the GUI
if ( nargin < 2 )
    guiMode = 0;
else
    guiMode = 1;
end

% Call the exportInfo method of the chebguiExporterEIG class, which takes care
% of extracting most information we need from GUIFILE.
expInfo = chebguiExporterEIG.exportInfo(guifile);

% Extract the needed fields from the EXPINFO struct.
dom = str2num(expInfo.dom);
allVarNames = expInfo.allVarNames;
indVarName = expInfo.indVarName;
eigVarName = expInfo.lname;
sigma = expInfo.sigma;
K = str2double(expInfo.K);
rhsString = expInfo.rhsString;

% Are we working with a generalized eigenvalue problem?
generalized = expInfo.generalized;

% Convert the LHS string to proper anonymous function using eval
LHS = eval(expInfo.lhsString);

% Support for periodic boundary conditions
if ( expInfo.periodic )
    BC = 'periodic';
else
    BC = eval(expInfo.bcString);
end

% Create the LHS chebop, and try to linearize it (if that gives an error, the
% operator is probably linear).
N_LHS = chebop(LHS, dom, BC);
try
    A = linop(N_LHS);
catch ME
    MEID = ME.identifier;
    if ( guiMode && ~isempty(strfind(MEID, 'linop:nonlinear')) )
        errordlg('Operator is not linear.', 'Chebgui error', 'modal');
    else
        rethrow(ME)
    end
    varargout{1} = handles;
    return
end

% Create the RHS chebop if the RHSSTRING is not empty, and try to linearize it
% (if that gives an error, the operator is probably linear).
if ( ~isempty(rhsString) )
    RHS = eval(rhsString);
    N_RHS = chebop(RHS, dom);

    try
        B = linop(N_RHS);
    catch ME
        MEID = ME.identifier;
        if ( guiMode  && ~isempty(strfind(MEID, 'linop:nonlinear'))  )
            errordlg('Operator is not linear.', 'Chebgui error', 'modal');
        else
            rethrow(ME)
        end
        varargout{1} = handles;
        return
    end
    
    % Did we determine that we had to negate the LHS operator?
    if ( expInfo.flipSigns )
        A = -A;
    end
    
end

% Obtain a CHEBOPPREF object
options = cheboppref;

% Check whether the tolerance is too tight.
defaultTol = options.errTol;
tolInput = guifile.tol;
if ( isempty(tolInput) )
    tolNum = defaultTol;
else
    tolNum = str2double(tolInput);
end

% Need to obtain a CHEBFUNPREF object to check what the current tolerance is set
% at.
chebfunp = chebfunpref;
if ( tolNum < chebfunp.techPrefs.eps )
    warndlg('Tolerance specified is less than current chebfun epsilon', ...
        'Warning','modal');
    uiwait(gcf)
end

% Do we want to show grid?
options.grid = guifile.options.grid;

% What discretization do we want?
options.discretization = expInfo.discretization;

% If we have specified the 'periodic' option, we need to throw away the boundary
% conditions from the LINOP A before continuing (as they're imposed by
% construction):
if ( strcmp(expInfo.discretization, 'periodic') )
    A.constraint = []; 
end

% Change various GUI components (only need to bother with in GUI mode).
if ( guiMode )
    set(handles.fig_sol, 'Visible', 'On');
    set(handles.fig_norm, 'Visible', 'On');
end

% Compute the eigenvalues!
if ( generalized )
    if ( isempty(sigma) )
        [V, D] = eigs(A, B, K, options);
    else
        [V, D] = eigs(A, B, K, sigma, options);
    end
else
    if ( isempty(sigma) )
        [V, D] = eigs(A, K, options);
    else
        [V, D] = eigs(A, K, sigma, options);
    end
end

% Sort the eigenvalues.
[D, idx] = sort(diag(D));

% We expect V to be a CHEBMATRIX. Sort the columns so they correspond to the now
% sorted eigenvalues in D.
V = V(:,idx);

if ( ~guiMode )
    % If we're not in GUI mode, we can finish here.
    if ( nargout == 1 )
        varargout = diag(D);
    else
        varargout{1} = D;
        varargout{2} = V;
    end
    
else
    % Now do some more stuff specific to GUI
    
    % Store in handles latest CHEBOP, eigenfunctions, eigenmodes etc. (enables
    % exporting later on):
    handles.latest.type = 'eig';
    handles.latest.solution = D;
    handles.latest.solutionT = V;
    handles.latest.chebop = A;
    handles.latest.options = options;
    
    % Notify the GUI we have a solution available
    handles.hasSolution = 1;
    handles.varnames = allVarNames;
    handles.eigVarName = eigVarName;
    handles.indVarName = indVarName{1};
    
    % Plot the eigenmodes.
    chebguiController.plotEigenmodes(handles, 0, handles.fig_sol, ...
        handles.fig_norm);
    
    set(handles.iter_text, 'Visible', 'on');
    set(handles.iter_text, 'String', 'Eigenvalues');
    set(handles.iter_list, 'Visible', 'on');
    
    % Display eigenvalues to level of tolerance
    s = num2str(ceil(-log10(tolNum)));
    set(handles.iter_list, 'String', num2str(D, ['%' s '.' s 'f']));
    set(handles.iter_list, 'Value', 1:numel(D));
    
    % Return the handles as varargout.
    varargout{1} = handles;
    
end

end

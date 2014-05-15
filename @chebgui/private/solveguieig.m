function varargout = solveguieig(guifile, handles)
% SOLVEGUIEIG

% TODO:  Documentation.

% Copyright 2014 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/chebfun/ for Chebfun information.

% Create a domain and the linear function on that domain. We use xt for the
% linear function, later in the code we will be able to determine whether x
% or t is used for the linear function.
defaultTol = 1e-10; % TODO: This should be read from cheboppref
% defaultTol = max(cheboppref('restol'),cheboppref('deltol'));

% Handles will be an empty variable if we are solving without using the GUI
if ( nargin < 2 )
    guiMode = 0;
else
    guiMode = 1;
end
dom = str2num(guifile.domain);
x = chebfun(@(x) x, dom);

% Extract information from the GUI fields
deInput = guifile.DE;
bcInput = guifile.BC;

sigma = [];
if ( ~isempty(guifile.sigma) )
    sigma = guifile.sigma;
    numSigma = str2num(sigma); 
    if ( ~isempty(numSigma) )
        sigma = numSigma;
    end
end
K = 6;
if ( isfield(guifile.options, 'numeigs') && ~isempty(guifile.options.numeigs) )
    K = str2double(guifile.options.numeigs);
end

% Wrap all input strings in a cell (if they're not a cell already)
if ( isa(deInput,'char') )
    deInput = cellstr(deInput);
end

if ( isa(bcInput,'char') )
    bcInput = cellstr(bcInput);
end

[allStrings allVarString indVarName ignored ignored eigVarName allVarNames] ...
    = setupFields(guifile, deInput, 'DE');
handles.varnames = allVarNames;

% If indVarName is empty, use the default value
if ( isempty(indVarName{1}) )
    indVarName{1} = 'x';
end
% Replace 'DUMMYSPACE' by the correct independent variable name
allStrings = strrep(allStrings, 'DUMMYSPACE', indVarName{1});
% Pretty print feval statements
allStrings = prettyprintfevalstring(allStrings, allVarNames);

% If allStrings return a cell, we have both a LHS and a RHS string. Else,
% we only have a LHS string, so we need to create the LHS linop manually.
if ( iscell(allStrings) )
    lhsString = allStrings{1};
    rhsString = allStrings{2};
else
    lhsString = allStrings;
    rhsString = '';
end
% Convert the strings to proper anon. function using eval
LHS = eval(lhsString);

% Support for periodic boundary conditions
BC = [];
if ( ~isempty(bcInput{1}) )
    if ( strcmpi(bcInput{1}, 'periodic') )
        bcInput{1} = [];
        BC = 'periodic';
    else
        bcString = setupFields(guifile, bcInput, 'BCnew', allVarString);
        bcString = strrep(bcString, 'DUMMYSPACE', indVarName{1});
        BC = eval(bcString);
    end
end

% Variable which determines whether it's a generalized problem. If
% rhsString is empty, we can be sure it's not a generalized problem.
generalized = 1;

% Create the chebops, and try to linearise them.
% We will always have a string for the LHS, if the one for RHS is empty, we
% know we have a non-generalised problem.
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
    
    % Check whether we are working with generalized
    % problems or not by comparing B with the identity operator on the domain.
    I = operatorBlock.eye(dom);
    % Set a discretization size for comparing operators
    discDim = 10;    
    % Obtain a discretisation of the operator B
    Bdisc = matrix(B, discDim);
    % Obtain a discretization of the identity operator on the domain
    Idisc = matrix(linop(I), discDim);
    % In case of systems, B will have a block structure. Need to get the
    % corresponding block identity operator by tiling the identity operator.
    % Start by wrapping in a cell and call repmat
    Idisc = repmat({Idisc}, [1, size(B, 1)]);
    % Expand to a block diagonal matrix
    Idisc = blkdiag(Idisc{:});
    
    % Compare the discretizations to see whether they are the same
    if ( size(Bdisc, 1) ~= size(Bdisc, 2) )
        generalized = 1;
    else

        opDifference = Bdisc - Idisc;
        opSum = Bdisc + Idisc;
        if ( isempty(nonzeros(opDifference)) )
            generalized = 0;
        end

        if ( isempty(nonzeros(opSum)) )
            generalized = 0;
            A = -A;
        end
    end
else
    generalized = 0;
end

tolInput = guifile.tol;
if ( isempty(tolInput) )
    tolNum = defaultTol;
else
    tolNum = str2double(tolInput);
end

chebfunp = chebfunpref;

if ( tolNum < chebfunp.techPrefs.eps )
    warndlg('Tolerance specified is less than current chebfun epsilon', ...
        'Warning','modal');
    uiwait(gcf)
end

options = cheboppref;
% Do we want to show grid?
options.grid = guifile.options.grid;

% Change various GUI components (only need to bother with in GUI mode).
if ( guiMode )

    % TODO:  Delete if no longer necessary.
%     set(handles.iter_list,'String','');
%     set(handles.iter_text,'Visible','On');
%     set(handles.iter_list,'Visible','On');

    set(handles.fig_sol, 'Visible', 'On');
    set(handles.fig_norm, 'Visible', 'On');
end

% Compute the eigenvalues.
if ( generalized )
    if ( isempty(sigma) )
        [V, D] = eigs(A, B, K);
    else
        [V, D] = eigs(A, B, K, sigma);
    end
else
    if ( isempty(sigma) )
        [V, D] = eigs(A, K);
    else
        [V, D] = eigs(A, K, sigma);
    end
end
[D, idx] = sort(diag(D));

% TODO: What happens in case of coupled systems.
% For now, assume scalar problems

% We expect V to be a CHEBMATRIX. Sort the columns so they correspond to the now
% sorted eigenvalues in D.
for k = 1:size(V, 1)
    V(k) = extractColumns(V{k}, idx);
end

% If we're not in GUI mode, we can finish here.
if ( ~guiMode )
    if ( nargout == 1 )
        varargout = diag(D);
    else
        varargout{1} = D;
        varargout{2} = V;
    end
    return
end

% Now do some more stuff specific to GUI
if ( guiMode )
    % Store in handles latest chebop, solution, vector of norm of updates etc.
    % (enables exporting later on)
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
    
    ploteigenmodes(handles.guifile, handles, 0, handles.fig_sol, ...
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

function str = prettyprintfevalstring(str, varnames)

for k = 1:numel(varnames)
    oldstr = ['feval(' varnames{k} ','];
    newstr = [varnames{k} '('];
    str = strrep(str, oldstr, newstr);
    oldstr = [varnames{k} '(''end'''];
    newstr = [varnames{k} '(end'];
    str = strrep(str, oldstr, newstr);
    oldstr = [varnames{k} '(''right'''];
    newstr = [varnames{k} '(end'];
    str = strrep(str, oldstr, newstr);
    oldstr = [varnames{k} '(''start'''];
    newstr = [varnames{k} '(' varnames{k} '.ends(1)'];
    str = strrep(str, oldstr, newstr);
    oldstr = [varnames{k} '(''left'''];
    newstr = [varnames{k} '(' varnames{k} '.ends(1)'];
    str = strrep(str, oldstr, newstr);
end

end

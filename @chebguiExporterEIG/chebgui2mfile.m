function chebgui2mfile(exporter, guifile, fid)
%EXPORTBVP2MFILE    Export a EIG problem from CHEBGUI to a .m file.
%
%   See also: chebgui/export.
% TODO:  Documentation.

% Copyright 2014 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/chebfun/ for Chebfun information.

% Extract information from the GUI fields
dom = guifile.domain;
deInput = guifile.DE;
bcInput = guifile.BC;

% Sigmas and num eigs
sigma = guifile.sigma;

% Convert to a proper string if needed
if ( ~isempty(sigma) && isempty(str2double(sigma)) )
    sigma = sprintf('''%s''',sigma);
end

K = guifile.options.numeigs;
if ( isempty(K) )
    K = '6';
end

% Wrap all input strings in a cell (if they're not a cell already)
if ( isa(deInput,'char') )
    deInput = cellstr(deInput);
end

if ( isa(bcInput,'char') )
    bcInput = cellstr(bcInput);
end

[allStrings allVarString indVarName pdeVarName pdeflag ignored allVarNames] ...
    = setupFields(guifile, deInput, 'DE');

% If indVarName is empty, use the default value
if ( isempty(indVarName{1}) )
    indVarName{1} = 'x';
end

% Replace 'DUMMYSPACE' by the correct independent variable name
allStrings = strrep(allStrings, 'DUMMYSPACE', indVarName{1});
% Pretty print feval statements
allStrings = exporter.prettyPrintFevalString(allStrings, allVarNames);

% If allStrings return a cell, we have both a LHS and a RHS string. Else,
% we only have a LHS string, so we need to create the LHS linop manually.
if ( iscell(allStrings) )
    lhsString = allStrings{1};
    rhsString = allStrings{2};
else
    lhsString = allStrings;
    rhsString = '';
end

% Assign x or t as the linear function on the domain
d = str2num(dom);
xt = chebfun('x', d);
eval([indVarName{1}, '=xt;']);

% Convert the strings to proper anon. function using eval
LHS  = eval(lhsString);

% Support for periodic boundary conditions
periodic = false;
BC = [];
if ( ~isempty(bcInput{1}) )
    if ( strcmpi(bcInput{1},'periodic') )
        bcInput{1} = [];
        periodic = true;
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
N_LHS = chebop(LHS, d, BC);

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
% Check for a generalised problem
if ( ~isempty(rhsString) )
    RHS  = eval(rhsString);
    N_RHS = chebop(RHS, d);

    try
        B = linop(N_RHS);
    catch ME
        MEID = ME.identifier;
        if ( guiMode  && ~isempty(strfind(MEID, 'linop:nonlinear')) )
            errordlg('Operator is not linear.', 'Chebgui error', 'modal');
        else
            rethrow(ME)
        end
        varargout{1} = handles;
        return
    end
    
    % Check whether we are working with generalized
    % problems or not by comparing B with the identity operator on the domain.
    I = operatorBlock.eye(str2num(dom));
    % Set a discretization size for comparing operators
    discDim = 10;    
    % Obtain a discretisation of the operator B
    Bdisc = matrix(B, discDim);
    % Obtain a discretization of the identity operator on the domain
    Idisc = matrix(linop(I), discDim);
    % In case of systems, B will have a block structure. Need to get the
    % corresponding block identity operator by tiling the identity operator.
    % Start by wrapping in a cell and call repmat
    Idisc = repmat({Idisc}, [1, size(B,1)]);
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

% Find the eigenvalue name
mask = strcmp(deInput{1}, {'lambda', 'lam', 'l'});

if ( mask(1) )
    lname = 'lambda'; 
elseif ( mask(2) )
    lname = 'lam'; 
elseif ( mask(3) )
    lname = 'l'; 
else
    lname = 'lambda';
end

% Print to EIG problem
fprintf(fid, '%% Solving\n');
if ( (numel(deInput) == 1) && ~any(deInput{1} == '=') )
            fprintf(fid, '%%   %s = %s*%s\n', deInput{1}, lname, allVarString);
else
    for k = 1:numel(deInput)
        fprintf(fid, '%%   %s,\n', deInput{k});
    end
end
fprintf(fid, '%% for %s in %s', indVarName{1}, dom);
if ( ~isempty(bcInput{1}) )
    fprintf(fid,',  subject to\n%%');
    for k = 1:numel(bcInput)
        fprintf(fid, '   %s', bcInput{k});
        if ( (k ~= numel(bcInput)) && (numel(bcInput) > 1) )
             fprintf(fid, ',\n%%');
        end
    end
    fprintf(fid, '.\n');
elseif ( periodic )
    fprintf(fid, ', subject to periodic boundary conditions.\n\n');
else
    fprintf(fid, '.\n');
end

fprintf(fid, '\n%%%% Define the domain we''re working on.\n');
fprintf(fid, 'dom = %s;\n', dom);
if ( ~generalized )
    fprintf(fid, ['\n%% Assign the equation to a chebop N such that' ...
        ' N(u) = %s*u.\n'], lname);
    fprintf(fid, 'N = chebop(%s, dom);\n', lhsString);
else
    fprintf(fid, ['\n%% Assign the equation to two chebops N and B such' ...
        'that N(u) = %s*B(u).\n'], lname);
    fprintf(fid, 'N = chebop(%s, dom);\n', lhsString);
    fprintf(fid, 'B = chebop(%s, dom);\n', rhsString);
end

% Make assignments for BCs.
fprintf(fid, '\n%% Assign boundary conditions to the chebop.\n');
if ( ~isempty(bcInput{1}) )
    bcString = exporter.prettyPrintFevalString(bcString, allVarNames);
    fprintf(fid, 'N.bc = %s;\n', bcString);
end
if ( periodic )
    fprintf(fid, 'N.bc = ''periodic'';\n');
end

fprintf(fid, '\n%% Number of eigenvalue and eigenmodes to compute.\n');
fprintf(fid, 'k = %s;\n', K);

fprintf(fid, '\n%% Number of eigenvalue and eigenmodes to compute.\n');
fprintf(fid, 'k = %s;\n', K);

fprintf(fid, '\n%%%% Solve the eigenvalue problem.\n');
if ( ~generalized )
    if ( ~isempty(sigma) )
        fprintf(fid, '[V, D] = eigs(N, k, %s);\n', sigma);
    else
        fprintf(fid, '[V, D] = eigs(N, k);\n');
    end
else
    if ( ~isempty(sigma) )
        fprintf(fid, '[V, D] = eigs(N, B, k, %s);\n', sigma);
    else
        fprintf(fid, '[V, D] = eigs(N, B, k);\n');
    end
end
fprintf(fid, '\n%%%% Plot the eigenvalues.\n');
fprintf(fid, 'D = diag(D);\n');
fprintf(fid, 'figure\n');
fprintf(fid, 'plot(real(D), imag(D), ''.'', ''markersize'', 25)\n');
fprintf(fid, 'title(''Eigenvalues''); xlabel(''real''); ylabel(''imag'');\n');

if ( ischar(allVarNames) || (numel(allVarNames) == 1) )
    fprintf(fid, '\n%% Plot the eigenmodes.\n');
    fprintf(fid, 'figure\n');
    fprintf(fid, 'plot(real(V), ''linewidth'', 2);\n');
    fprintf(fid, 'title(''Eigenmodes''); xlabel(''%s''); ylabel(''%s'');\n', ...
        indVarName{1}, allVarString);
end

end
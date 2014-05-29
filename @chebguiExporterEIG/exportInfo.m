function expInfo = exportInfo(e, guifile)

% Extract information from the CHEBGUI fields
dom = guifile.domain;
deInput = guifile.DE;
bcInput = guifile.BC;


% Wrap all input strings in a cell (if they're not a cell already)
if ( isa(deInput,'char') )
    deInput = cellstr(deInput);
end

if ( isa(bcInput,'char') )
    bcInput = cellstr(bcInput);
end

% The sigma option for specifying what eigenvalues we seek:
sigma = guifile.sigma;

% Convert to a proper string if needed
if ( ~isempty(sigma) && isempty(str2double(sigma)) )
    sigma = sprintf('''%s''',sigma);
end

% Number of eigenvalues we want (default value: 6):
K = guifile.options.numeigs;
if ( isempty(K) )
    K = '6';
end

[allStrings, allVarString, indVarName, dummy, dummy, dummy, allVarNames] ...
    = setupFields(guifile, deInput, 'DE');

% If indVarName is empty, use the default value
if ( isempty(indVarName{1}) )
    indVarName{1} = 'x';
end

% Replace 'DUMMYSPACE' by the correct independent variable name
allStrings = strrep(allStrings, 'DUMMYSPACE', indVarName{1});
% Pretty print feval statements
allStrings = e.prettyPrintFevalString(allStrings, allVarNames);

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
LHS = eval(lhsString);

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
% Check for a generalised problem.
% TODO: This should be a joint method with solveGUIeig
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

%% Fill up the expInfo struct
expInfo.dom = dom;
expInfo.deInput = deInput;
expInfo.bcInput = bcInput;

% Eigenvalue problem specific information
expInfo.K = K;
expInfo.sigma = sigma;
expInfo.generalized = generalized;
expInfo.lname = lname;


expInfo.lhsString = lhsString;
expInfo.allStrings = allStrings;
expInfo.allVarString = allVarString;
expInfo.allVarNames = allVarNames;
expInfo.indVarName = indVarName;
expInfo.periodic = periodic;
expInfo.bcString = bcString;

% Information related to options set-up
expInfo.discretization = guifile.options.discretization;
end
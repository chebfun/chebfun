function expInfo = exportInfo(guifile)
%EXPORTINFO   Extract useful info from a CHEBGUI object for exporting
%
% Calling sequence
%
%   EXPINFO = EXPORTINFO(GUIFILE)
%
% where
%
%   GUIFILE:    A CHEBGUI object
%   EXPINFO:    A struct, containing fields with information for exporting to an
%               .m-file.

% Copyright 2015 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

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

% Find the eigenvalue parameter name:
mask = cellfun(@strfind, repmat(deInput(1), 1, 3), {'lambda', 'lam', 'l'}, ...
    'UniformOutput', false);
if ( ~isempty(mask{1}) )
    lname = 'lambda'; 
elseif ( ~isempty(mask{2}) )
    lname = 'lam'; 
elseif ( ~isempty(mask{3}) )
    lname = 'l'; 
else
    lname = '';
end

% Ensure that l, lam or lambda appears in the problem!
assert(~isempty(lname), 'CHEBFUN:CHEBGUIEXPORTEREIG:exportInfo:lname', ...
    ['Variable representing eigenvalue parameter not found in Chebgui input. ' ...
    'Please ensure ''lambda'', ''lam'' or ''l'' appears in input.'])

% Obtain strings for setting up the problem:
[allStrings, allVarString, indVarName, dummy1, pdeflag, dummy3, allVarNames] ...
    = setupFields(guifile, deInput, 'DE');

% If indVarName is empty, use the default value
if ( isempty(indVarName{1}) )
    indVarName{1} = 'x';
end

% Replace 'DUMMYSPACE' by the correct independent variable name
allStrings = strrep(allStrings, 'DUMMYSPACE', indVarName{1});
% Pretty print feval statements
allStrings = chebguiExporter.prettyPrintFevalString(allStrings, allVarNames);

% If allStrings return a cell, we have both a LHS and a RHS string. Else,
% we only have a LHS string, so we need to create the LHS linop manually.
if ( iscell(allStrings) )
    lhsString = allStrings{1};
    rhsString = allStrings{2};
else
    lhsString = allStrings;
    rhsString = '';
end

% Assign x or t as the linear function on the domain. This is required so that
% we can check whether the problem is a generalized eigenvalue problem or not:
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
        bcString = [];
    else
        bcString = setupFields(guifile, bcInput, 'BCnew', allVarString);
        bcString = strrep(bcString, 'DUMMYSPACE', indVarName{1});
        BC = eval(bcString);
    end
end

% What discretization option do we want?
discretization = chebguiExporter.discOption(periodic, dom, ...
    guifile.options.discretization);

% We need a cheboppref to be able to call matrix() below to determine whether
% the problem is generalized.
pref = cheboppref();
pref.discretization = discretization;

% Variable which determines whether it's a generalized problem. If
% rhsString is empty, we can be sure it's not a generalized problem.
generalized = 1;

% Check for a generalised problem.
if ( ~isempty(rhsString) )
    RHS  = eval(rhsString);
    N_RHS = chebop(RHS, d);

    % Linearize the RHS operator, so that it can be compared with the identity
    % operator.
    linCheck = false;
    paramReshape = false;
    B = linearize(N_RHS, [], [], linCheck, paramReshape);
    
    % Set a discretization size for comparing operators
    discDim = repmat(10, 1, length(d) - 1);    
    % Obtain a discretisation of the operator B
    Bdisc = matrix(B, discDim, pref);
    % Obtain a discretization of the chebop specified via chebop(@(u) u). Notice
    % that this is not identical to the discretization of identity operator on
    % the domain, due to the fact that chebop/linearize now tries to
    % automatically detect parameters appearing in the problem!
    Idisc = ones(10*(length(d) - 1), 1);
    % In case of systems, B will have a block structure. Need to get the
    % corresponding block identity operator by tiling the identity operator.
    % Start by wrapping in a cell and call repmat
    Idisc = repmat({Idisc}, [1, size(B,1)]);
    % Expand to a block diagonal matrix
    Idisc = blkdiag(Idisc{:});
    
    % Compare the discretizations to see whether they are the same. If Bdisc is
    % not tall and skinny, it can't have been the identity operator! In other
    % words, if Bdisc has more columns than number of variables that appear in
    % the problem, it will not have been the identity operator.
    if ( size(Bdisc, 2) ~= size(B, 2) )
        generalized = 1;
    else
        
        opDifference = Bdisc - Idisc;
        opSum = Bdisc + Idisc;
        tol = 10*(length(d) - 1)*size(B, 2)*eps;
        % Check whether B matches identity (allow for numerical roundoff errors)
        if ( norm(opDifference) < tol )
            generalized = 0;
        end

        if ( norm(opSum) < tol )
            generalized = 0;
        end
    end
else
    generalized = 0;
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

% And then some...
expInfo.lhsString = lhsString;
expInfo.rhsString = rhsString;
expInfo.allStrings = allStrings;
expInfo.allVarString = allVarString;
expInfo.allVarNames = allVarNames;
expInfo.indVarName = indVarName;
expInfo.periodic = periodic;
expInfo.bcString = bcString;

% Information related to options set-up
expInfo.discretization = discretization;

end

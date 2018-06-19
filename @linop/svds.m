function varargout = svds(A, varargin)
%SVDS   Find selected singularvalues and singularfunctions of a linear operator.
%   Important (1): While you can construct a LINOP and apply this method, the
%   recommended procedure is to use CHEBOP/SVDS instead.
%   Important (2): A CHEBOPPREF object PREFS has to be passed. When this method
%   is called via CHEBOP/SVDS, PREFS is inherited from the CHEBOP level.
%   Important (3): This method depends on LINOP/ADJOINT and can only be used
%   for LINOPS representing differential operators. Strictly multiplication,
%   integral or other exotic operators are not supported.
%
%   S = SVDS(A, PREFS) returns a vector of 6 singularvalues of the linop A.
%   SVDS will attempt to return the singularvalues corresponding to the most
%   easily resolved singularfunctions. (This is unlike the built-in SVDS, which
%   returns the largest singularvalues by default.)
%
%   [U, S, V] = SVDS(A, PREFS) returns a diagonal 6x6 matrix S of A's most
%   easily resolved singularvalues, and their corresponding left and right
%   singularfunctions in the chebmatrices U and V respectively, where V{i}(:,j)
%   is the jth singularfunction in variable i of the system.
%
%   [...] = SVDS(A, K, PREFS) find the K most easily resolved 
%   singularvalues.
%
%   This version of SVDS does not use iterative methods as in the built-in
%   SVDS for sparse matrices. Instead, it uses the built-in EIG on dense
%   matrices of increasing size, stopping when the targeted singularfunctions
%   appear to have converged, as determined by the chebfun constructor.
%
%   EXAMPLE: First derivative operator
%
%   d = [0 pi];
%   A = linop( operatorBlock.diff(d) );
%   prefs = cheboppref();
%   prefs.discretization = @chebcolloc2;
%   [U,S,V] = svds(A, 10, prefs);
%   format long, sqrt(-diag(D))  % integers, to 14 digits
%
% See also CHEBOP/SVDS, CHEBOP/ADJOINT, LINOP/EIGS and LINOP/ADJIONT.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Parsing inputs.
k = [];       % will be made default value below
prefs = [];
bcType = [];
gotk = false; % until we detect a value of k in inputs
for j = 1:nargin-1
    item = varargin{j};
    if ( isa(item,'cheboppref') )
        prefs = item;
    elseif ( ~gotk && isnumeric(item) && (item > 0) && (item == round(item) ) )
        k = item;
        gotk = true;
    elseif ( isa(item,'char') )
        bcType = item;
    else
        error('Could not parse argument number %i.',j+1)
    end
end

% Check that we recieved a prefs object
if ( isempty(prefs) )
    error('CHEBFUN:LINOP:svds:prefs', ...
        'A preference object is required.');
end

% Check that we recieved a bcType
if ( isempty(bcType) )
    error('CHEBFUN:LINOP:svds:bcType', ...
        'A bcType is required.');
end

% Check for unbounded domains:
if ( ~all(isfinite(A.domain)) )
    error('CHEBFUN:LINOP:eigs:infDom', ...
        'Unbounded domains are not supported.');
end

% Assign default to k if needed.
if ( isempty(k) || isnan(k) )
    k = 6;
end

% Construct adjoint
Astar = linopAdjoint(A,bcType);

% Initialize superA to 0
[m, n] = size(A);
dom = A.domain;
nm = n + m;
superA = linop(mat2cell(zeros(nm),ones(1,nm),ones(1,nm)));
z = chebfun(0,dom);
for ii = 1:nm
    for jj = 1:nm
        superA.blocks{ii,jj} = operatorBlock.mult(z,z.domain);
    end
end
superA.domain = A.domain;

% Fill the appropriate blocks
for ii = 1:m
    for jj = 1:n
        superA.blocks{n+ii,jj} = A.blocks{ii,jj};
    end
end
for ii = 1:n
    for jj = 1:m
        superA.blocks{ii,n+jj} = Astar.blocks{ii,jj};
    end
end

% Get constraints
C = A.constraint;
Cstar = Astar.constraint;

% Get funs
Cfuns = C.functional; nc = size(Cfuns,1);
Csfuns = Cstar.functional; ncs = size(Csfuns,1);

% Make superC
dor = max(max(A.diffOrder));
superC.functional = [];
z = functionalBlock.zero(dom);
for ii = 1:nc
    row = [];
    for jj = 1:n
        row = [row,Cfuns{ii,jj}];
    end
    for jj = 1:m
        row = [row,z];
    end
    superC.functional = [superC.functional;row];
end
for ii = 1:ncs
    row = [];
    for jj = 1:n
        row = [row,z];
    end
    for jj = 1:m
        row = [row,Csfuns{ii,jj}];
    end
    superC.functional = [superC.functional;row];
end
superC.values = zeros(dor*nm,1);

% Finish superA
superA.constraint = superC;

% Count number of null vectors for A and Astar
nulA = dor*n-nc;

% Set number of singular values. If the null space is empty then there are
% exactly 2 copies +/- of k singular values. If the null space is non-empty then
% those zero singular values won't be repeated. We compute two extra just to
% make sure got didn't lose one to a sign change.
nsvals = 2 + 2*k-abs(nulA);

% Call linop/eigs using sigma = 0 since we are assuming L is an unbounded
% differential operator
% warning('off','all') % turn warnings off
if ( strcmp(bcType,'periodic') )
    prefs.discretization = @trigcolloc;
else
    prefs.discretization = @chebcolloc2;
end
[ Q, D ] = eigs( superA, nsvals, 0, prefs, 'rayleigh' );
% warning('on','all') % turn warnings back on

% make sure singular values are real
if ( any(imag(diag(D)) ~= 0) )
    error('CHEBFUN:LINOP:svds:real', ...
        'Computed singular values are not strictly real.');
end

% Set tiny values of D to zero
D = diag(D);
D(abs(D) < prefs.bvpTol*max(abs(D))) = 0;

% Sort by first inverting the singular values. This works because for
% differential operators the smoothest singular functions always correspond to
% the tiniest singular values.
[D,id] = sort(1./D,'descend');
Q = Q(:,id);

% Flip, reorder and discard unwanted functions.
S = diag(1./D(k:-1:1));
Q = Q(:,k:-1:1);

% Rescale singular vectors
V = Q(1:n,:); nrmV = sqrt(diag(V'*V)); 
nrmV( nrmV < prefs.bvpTol ) = 1; 
V = V*diag(1./nrmV);
U = Q(n+1:end,:); 
nrmU = sqrt(diag(U'*U));
nrmU( nrmU < prefs.bvpTol ) = 1; 
nrmU = 1./nrmU;
U = U*diag(nrmU);

% Set output
if ( nargout <= 1 )
    varargout = { diag(S) };
elseif ( nargout == 3 )
    varargout = { U, S, V };
else
    error('CHEBFUN:LINOP:svds:outputs', 'SVDS requires one or three outputs.');
end

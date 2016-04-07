function varargout = eigs(A, varargin)
%EIGS    Eigenvalues and eigenfunctions of a linear operator.
%   Important (1): While you can construct a LINOP and apply this method, the
%   recommended procedure is to use CHEBOP/EIGS instead.
%   Important (2): A CHEBOPPREF object PREFS has to be passed. When this method
%   is called via CHEBOP/EIGS, PREFS is inherited from the CHEBOP level.
%
%   D = EIGS(A, PREFS) returns a vector of 6 eigenvalues of the linop A. EIGS 
%   will attempt to return the eigenvalues corresponding to the most easily
%   resolved eigenfunctions. (This is unlike the built-in EIGS, which
%   returns the largest eigenvalues by default.)
%
%   [V, D] = EIGS(A, PREFS) returns a diagonal 6x6 matrix D of A's most easily
%   resolved eigenvalues, and their corresponding eigenfunctions in the
%   chebmatrix V, where V{i}(:,j) is the jth eigenfunction in variable i of
%   the system.
%
%   [...] = EIGS(A, B, PREFS) solves the generalized eigenproblem A*V = B*V*D,
%   where B is another linop.
%
%   EIGS(A, K, PREFS) and EIGS(A, B, K, PREFS) find the K most easily resolved 
%   eigenvalues.
%
%   EIGS(A, K, SIGMA, PREFS) and EIGS(A, B, K, SIGMA, PREFS) find K eigenvalues. 
%   If SIGMA is a scalar, the eigenvalues found are the ones closest to SIGMA. 
%   Other selection possibilities for SIGMA are:
%
%      'LM' (or Inf) and 'SM' for largest and smallest magnitude
%      'LR' and 'SR' for largest and smallest real part
%      'LI' and 'SI' for largest and smallest imaginary part
%
%   SIGMA must be chosen appropriately for the given operator. For example,
%   'LM' for an unbounded operator will fail to converge.
%
%   [...] = EIGS(A, ..., 'rayleigh') performs one step of Rayleigh quotient
%   iteration on the computed eigenpairs in an attempt to improve accuracy.
%
%   This version of EIGS does not use iterative methods as in the built-in
%   EIGS for sparse matrices. Instead, it uses the built-in EIG on dense
%   matrices of increasing size, stopping when the targeted eigenfunctions
%   appear to have converged, as determined by the chebfun constructor.
%
%   EXAMPLE: Simple harmonic oscillator
%
%   d = [0 pi];
%   A = linop( operatorBlock.diff(d, 2) );
%   E = functionalBlock.eval(d);
%   A = addBC(A, E(0), 0);
%   A = addBC(A, E(pi), 0);
%   prefs = cheboppref();
%   prefs.discretization = @chebcolloc2;
%   [V,D] = eigs(A, 10, prefs);
%   format long, sqrt(-diag(D))  % integers, to 14 digits
%
% See also CHEBOPPREF, CHEBOP.EIGS.

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Parsing inputs.
B = [];       % no generalized operator
k = [];       % will be made default value below
sigma = [];   % default 'auto' mode
prefs = [];
rayleigh = false; % do not perform rayleigh quotient iteration by default
gotk = false; % until we detect a value of k in inputs
for j = 1:nargin-1
    item = varargin{j};
    if ( isa(item, 'linop') )
        % Generalized operator term
        B = item;
    elseif ( isa(item,'cheboppref') )
        prefs = item;
    elseif ( ~gotk && isnumeric(item) && (item > 0) && (item == round(item) ) )
        % k should be given before sigma (which might also be integer)
        k = item;
        gotk = true;
    elseif ( strcmpi(item,'rayleigh') )
        rayleigh = true;
    elseif ( ischar(item) || isnumeric(item) )
        sigma = item;
    else
        error('Could not parse argument number %i.',j+1)
    end
end

% Check for unbounded domains:
if ( ~all(isfinite(A.domain)) )
    error('CHEBFUN:LINOP:eigs:infDom', ...
        'Unbounded domains are not supported.');
end

%#ok<*ASGLU> % Prevent MLINT warnings for unused variables, which are used in 
             % many places in this code to avoid the [~, arg2] = ... syntax.
             
% Discretization type:
discType = prefs.discretization;

% Make sure we have a valid discretization preference at this level.
assert(~ischar(discType), 'CHEBFUN:LINOP:expm:discretization', ...
    'pref.discretization must be a function handle, not a string.');

% Assign default to k if needed.
if ( isempty(k) || isnan(k) )
    k = 6;
end

% Check for square operator. (This is not strict enough, technically.)
m = size(A, 2);
if ( m ~= size(A, 1) )
    error('CHEBFUN:LINOP:eigs:notSquare','Block size must be square.')
end


% If there is a generalized eigenproblem, the domains must be merged before
% deriving the continuity equations for A:
if ( ~isempty(B) )
    % Merge the domains of A and B:
    dom = domain.merge(A.domain, B.domain);
    A.domain = dom;
    B.domain = dom;
end

% Set up the discretization of A:
if ( isa(discType, 'function_handle') )
    % Create a discretization object
    discA = discType(A);

    % Set the allowed discretisation lengths:
    dimVals = discA.dimensionValues(prefs);

    % Update the discretiztion dimension on unhappy pieces:
    discA.dimension = repmat(dimVals(1), 1, numel(discA.domain)-1);
else
    % A discretization is given:
    discA = discType;

    % Initialise dimVals;
    dimVals = max(discA.dimension);
end

% Construct a discretization for B:
if ( ~isempty(B) )

    constructor = str2func( class(discA) );   % constructor handle.
    discB = constructor(B);
    
    % We can ignore constraints and continuity--enforced on the left side.
    if ( ~isempty(discB.source.constraint) )
        discB.source.constraint = [];
        warning('CHEBFUN:LINOP:eigs:constraints', ...
                'Constraints on B are ignored.')
    end
    if ( ~isempty(discB.source.continuity) )
        discB.source.continuity = [];
        warning('CHEBFUN:LINOP:eigs:continuity', ...
                'Continuity conditions on B are ignored.')
    end       
    
    % Merge the two discretizations:
    [discA, discB] = merge(discA, discB);
    
else
    discB = [];
end

if ( isempty(A.continuity) )
     % Apply continuity conditions:
     discA.source = deriveContinuity(discA.source);
end

% 'SM' is equivalent to eigenvalues nearest zero.
if ( strcmpi(sigma, 'SM') )
    sigma = 0;
end

% Information required for finding the eigenvalues and functions.
numInts = discA.numIntervals;
numVars = size(A,2);

% compute adaptive grid sizes
minDE = max(floor(log2(prefs.minDimension)),ceil(log2(k)));
maxDE = ceil(log2(prefs.maxDimension));
Dims = 2.^(minDE:maxDE);
Dims(1) = max(Dims(1),prefs.minDimension);
Dims(end) = min(Dims(end),prefs.maxDimension);

% Loop through grids until converged.
for dim = Dims
     
    % set the dimension
    discA.dimension = dim*ones(1, numInts);

    % Discretize the LHS operator (incl. constraints/continuity):
    [PA, P, C, ignored, PS] = matrix(discA);
    
    % Discretize the RHS operator, or use identity.
    if ( ~isempty(discB) )
        % TODO: This is untidy. Can we make a method to do this? NH Apr 2014.
        discB.dimension = discA.dimension;
        PB = matrix(discB);
        % Project RHS matrix and prepend rows for the LHS constraints.
        PB = [ zeros(size(C)) ; PB ];
    else
        PB = [ zeros(size(C)) ; PS ];
    end    

    % call eig
    [V, D] = eig(full(PA), full(PB));

    % remove infinite eigenvalues
    [~,idx] = sort(abs(diag(D)),'descend');
    idx = idx(size(C)+1:end);
    D = D(idx,idx);
    V = V(:,idx);

    % Convert the discrete V values to CHEBMATRIX
    u = mat2fun(discA,P*V);
    u = vertcat(u{:});

    % simplify each chebfun
    for jj = 1:size(u,2)
        u(:,jj) = simplify(u(:,jj),prefs.bvpTol);
    end
    
    % compute sum of lengths of each piece in each column
    lens = cellfun(@length,{u{1:end,1:end}});
    lens = max(reshape(lens,size(u)),[],1);

    % sort eigenvalues by user prescribed flag
    if ( ~isempty(sigma) )

        % sort eigenvalues and eigenfunctions
        inds = nearest(diag(D),sigma);

    % sort according to smoothness
    else

        % sort converged functions by length
        inds = 1:length(lens);
        [~,idx] = sort(lens(inds),'ascend');
        inds = inds(idx);

    end

    % sort data according to inds
    inds = inds(1:k);
    lens = lens(inds)
    D = D(inds,inds);
    u = u(:,inds);
 
    % only proceed if at least k functions have converged
    if ( all(lens+5 < dim*numInts) )

        % use default matlab sort
        [lam,inds] = sort(diag(D));
        D = D(inds,inds);
        u = u(:,inds);

        % normalize
        for ii = 1:k
            u(:,ii) = u(:,ii)/norm(u(:,ii));
        end

        % do one step of Rayleigh quotient iteration to improve accuracy
        if ( rayleigh ) 
            [u, D] = rayleighQI(A,B,u,D,prefs);
        end

        % one output
        if ( nargout <= 1 )
            varargout = {lam};
        % two outputs
        elseif ( nargout == 2 )
            varargout = {u, D};
        % more than two outputs
        else
            error('CHEBFUN:linop:eigs','Maximum of two outputs.');
        end
  
        % return
        return

    end

    % throw error in max its hit
    if ( dim == prefs.maxDimension )
	error('CHEBFUN:linop:eigs',...
              'Maximum dimension reached without convergence.');
    end    

end

end
% END OF MAIN FUNCTION

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function idx = nearest(lam, sigma)
% Returns index vector that sorts eigenvalues by the given criterion.

if ( isempty(sigma) )
    idx = 1:length(lam);
elseif ( isnumeric(sigma) )
    if ( isinf(sigma) )
        [junk, idx] = sort(abs(lam), 'descend');
    else
        [junk, idx] = sort(abs(lam-sigma), 'ascend');
    end
else
    switch upper(sigma)
        case 'LR'
            [junk, idx] = sort(real(lam), 'descend');
        case 'SR'
            [junk, idx] = sort(real(lam), 'ascend');
        case 'LI'
            [junk, idx] = sort(imag(lam), 'descend');
        case 'SI'
            [junk, idx] = sort(imag(lam), 'ascend');
        case 'LM'
            [junk, idx] = sort(abs(lam), 'descend');
        case 'SM'
            [junk, idx] = sort(abs(lam), 'ascend');
        otherwise
            error('CHEBFUN:LINOP:eigs:sigma', 'Unidentified input ''sigma''.');
    end
end

% Delete infinite values. These can arise from rank deficiencies in the
% RHS matrix of the generalized eigenproblem.
idx( ~isfinite(lam(idx)) ) = [];

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [u, D] = rayleighQI(A, B, u, D, prefs)
% Function to perform one step of Rayleigh quotient iteration.

    % If B isempty set to identity
    if ( isempty(B) )
        B = 0*A;
        I = operatorBlock.eye(A.domain);
        for ii = 1:size(B,1)
            B.blocks{ii,ii} = I;
        end
    end
    B = linop(B);

    % Compute current Rayleigh quotients
    d = diag(D);

    % Loop through d
    for ii = 1:length(d)
        lam = d(ii);
        L = linop(A - lam*B);
        L.constraint = A.constraint;
        rhs = B*u(:,ii);
        v = linsolve(L, rhs, prefs);
        u(:,ii) = v/norm(v);
    end

    % Update D
    d = diag(u'*(A*u)) ./ diag(u'*(B*u));
    D = diag(d);

end

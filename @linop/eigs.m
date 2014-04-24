function varargout = eigs(L,varargin)
%EIGS    Eigenvalues and eigenfunctions of a linear operator.
%   Important: While you can construct a LINOP and apply this method, the
%   recommended procedure is to use CHEBOP.EIGS instead.
%
%   D = EIGS(A) returns a vector of 6 eigenvalues of the linop A. EIGS will
%   attempt to return the eigenvalues corresponding to the most easily
%   resolved eigenfunctions. (This is unlike the built-in EIGS, which
%   returns the largest eigenvalues by default.)
%
%   [V,D] = EIGS(A) returns a diagonal 6x6 matrix D of A's most easily
%   resolved eigenvalues, and their corresponding eigenfunctions in the
%   chebmatrix V, where V{i}(:,j) is the jth eigenfunction in variable i of
%   the system.
%
%   [...] = EIGS(A,B) solves the generalized eigenproblem A*V = B*V*D,
%   where B is another linop.
%
%   EIGS(A,K) and EIGS(A,B,K) find the K most easily resolved eigenvalues.
%
%   EIGS(A,K,SIGMA) and EIGS(A,B,K,SIGMA) find K eigenvalues. If SIGMA is
%   a scalar, the eigenvalues found are the ones closest to SIGMA. Other
%   selection possibilities for SIGMA are:
%
%      'LM' (or Inf) and 'SM' for largest and smallest magnitude
%      'LR' and 'SR' for largest and smallest real part
%      'LI' and 'SI' for largest and smallest imaginary part
%
%   SIGMA must be chosen appropriately for the given operator. For example,
%   'LM' for an unbounded operator will fail to converge.
%
%   EIGS(..., PREFS) accepts a CHEBOPPREF to control the behavior of
%   the algorithm. If empty, defaults are used.
%
%   This version of EIGS does not use iterative methods as in the built-in
%   EIGS for sparse matrices. Instead, it uses the built-in EIG on dense
%   matrices of increasing size, stopping when the targeted eigenfunctions
%   appear to have converged, as determined by the chebfun constructor.
%
%   EXAMPLE: Simple harmonic oscillator
%
%   d = [0 pi];
%   A = linop( operatorBlock.diff(d,2) );
%   E = functionalBlock.eval(d);
%   A = addBC(A,E(0),0);
%   A = addBC(A,E(pi),0);
%   [V,D] = eigs(A,10);
%   format long, sqrt(-diag(D))  % integers, to 14 digits
%
%   See also CHEBOPPREF, CHEBOP.EIGS.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org for Chebfun information.

% Parsing inputs.
M = [];       % no generalized operator
k = [];       % will be made default value below
sigma = [];   % default 'auto' mode
prefs = [];
gotk = false; % until we detect a value of k in inputs
for j = 1:nargin-1
    item = varargin{j};
    if ( isa(item, 'linop') )
        % Generalized operator term
        M = item;
    elseif ( isa(item,'cheboppref') )
        prefs = item;
    elseif ( ~gotk && isnumeric(item) && (item > 0) && (item == round(item) ) )
        % k should be given before sigma (which might also be integer)
        k = item;
        gotk = true;
    elseif ( ischar(item) || isnumeric(item) )
        sigma = item;
    else
        error('Could not parse argument number %i.',j+1)
    end
end

% Grab defaults if needed.
if ( isempty(prefs) )
    prefs = cheboppref;
end

% Discretization type.
discType = prefs.discretization;

% Assign default to k if needed.
if ( isempty(k) || isnan(k) )
    k = 6;
end

% Check for square operator. (This is not strict enough, technically.)
m = size(L, 2);
if ( m ~= size(L, 1) )
    error('LINOP:eigs:notsquare','Block size must be square.')
end

%% Set up the discretization:
if ( isa(discType, 'function_handle') )
    % Create a discretization object
    disc = discType(L);

    % Set the allowed discretisation lengths:
    dimVals = prefs.dimensionValues;

    % Update the discretiztion dimension on unhappy pieces:
    disc.dimension = repmat(dimVals(1), 1, numel(disc.domain)-1);
    dimVals(1) = [];
else
    % A discretization is given:
    disc = discType;

    % Initialise dimVals;
    dimVals = max(disc.dimension);
end

% Attach a domain to the discretization
dom = L.domain;
disc.domain = dom;

if ( isempty(L.continuity) )
     % Apply continuity conditions:
     disc.source = deriveContinuity(disc.source);
end

% If there is a generalized eigenproblem, the right-side operator needs to have
% its domain merged in and its own discretization.
discM = [];
if ( ~isempty(M) )
    dom = chebfun.mergeDomains(disc.domain,dom,M.domain);
    disc.domain = dom;   % update the discretization domain for L
    constructor = str2func( class(disc) );   % constructor handle
    discM = constructor(M,disc.dimension,disc.domain);
    % We can ignore constraints and continuity--enforced on the left side.
end


% 'SM' is equivalent to eigenvalues nearest zero.
if ( strcmpi(sigma, 'SM') )
    sigma = 0;
end

% Information required for finding the eigenvalues and functions.
numInts = disc.numIntervals;
isFun = isFunVariable(L);

% Automatic mode: find the best sigma by going where the convergence appears to
% be fastest.
if ( isempty(sigma) )
    % Try to determine where the 'most interesting' eigenvalue is.
    disc.dimension = 33*ones(1, numInts);
    [V1,D1] = getEigenvalues(disc, discM, 33, 0);
    disc.dimension(:) = 65;
    [V2,D2, P] = getEigenvalues(disc, discM, 33, 0);
    lam1 = diag(D1);
    lam2 = diag(D2);
    dif = bsxfun(@minus, lam1.', lam2);
    delta = min( abs(dif) );   % diffs from 33->65
    % These are the significantly large differences from 33->65.
    bigDel = (delta > 1e-12*norm(lam1,Inf));

    % Trim off things that are still changing a lot (relative to new size).
    lam1b = lam1;
    lam1b(bigDel) = 0;
    bigDel = logical((delta > 1e-3*norm(lam1b,Inf)) + bigDel);

    if ( all(bigDel) )
        % All values changed somewhat-- choose the one changing the least.
        [~,idx] = min(delta);
        sigma = lam1(idx);
    else
        % One by one, convert the eigenvectors to functions and check their cheb
        % expansion coefficients.
        U = partition(disc, P*V2);  % each cell is array valued, for one variable

        % Combine the different variable components into a single variable for
        % coefficient conversion.
        Z = 0;
        for j = ( find(isFun) )
            Z = Z + U{j};
        end

        % Convert the discrete Z values to CHEBFUN
        z = toFunction(disc, Z);

        % Obtain all coefficients to use below
        coeffs = get(z, 'coeffs');
        
        % Compute the 1-norm of the polynomial expansions, summing over smooth
        % pieces, for all columns.
        onenorm = 0;
        for j = 1:disc.numIntervals
            onenorm = onenorm + sum( abs( coeffs{j} ), 1 )';
        end
        
        [~, index] = min(onenorm);
        sigma = lam2(index);
    end
end

% Linear combination coefficients for convergence test. The convergence of the
% combination is the same as the worst constituent function. The nontrivial
% coefficents are to make accidental cancellations extremely unlikely.
coeff = 1./(2*(1:k)');

for dim = dimVals

    [V,D, P] = getEigenvalues(disc, discM, k, sigma);

    % Combine the eigenfunctions into a composite.
    v = V*coeff(1:size(V,2));

    % Convert the different components into cells
    u = partition(disc,P*v);

    % Test the happieness of the function pieces:
    [isDone, epsLevel] = testConvergence(disc, u(isFun));

    if ( all(isDone) )
        break
    else
        % Update the discretiztion dimension on unhappy pieces:
        disc.dimension(~isDone) = dim;
    end

end

% Detect finite rank operators.
if ( size(D,1) < k )
    if ( gotk )
        warning('CHEBFUN:linop:eigs:rank',...
            'Input has finite rank, only %d eigenvalues returned.', size(D,1));
    end
    k = size(D,1);
end

if ( nargout < 2 )  % Return the eigenvalues only
    varargout = { diag(D) };
else            % Unwrap the eigenvectors for output

    u = mat2fun(disc, V);

    % Find the norm in each eigenfunction (aggregated over variables).
    nrmsq = zeros(1,k);
    for j = 1:length(u)
        if ( isFun(j) )
            % Compress the representation.
            u{j} = simplify( u{j}, max(eps,epsLevel) );
        end
        nrmsq = nrmsq + sum( u{j}.*conj(u{j}), 1 );
    end

    % Normalize each eigenfunction.
    scale = diag( 1./sqrt(nrmsq') );
    for j = 1:length(u)
        u{j} = u{j}*scale;
    end

     varargout = { chebmatrix(u), D };
end

end
% END OF MAIN FUNCTION

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
function [V,D,P] = getEigenvalues(disc, discM, k, sigma)
% Formulate the discrete problem and solve for the eigenvalues

    % Discretize the LHS operator (incl. constraints/continuity):
    [PA, P, C, A] = matrix(disc);

    % Discretize the RHS operator, or use identity.
    if ( ~isempty(discM) )
        discM.dimension = disc.dimension;
        [~,~,~,B] = matrix(discM);
        % Project RHS matrix and prepend rows for the LHS constraints.
        PB = [ zeros(size(C)); P*B ];
    else
        PB = [ zeros(size(C)); P ];
    end

    % Compute eigenvalues.
    if ( length(PA) <= 2000 )
        [V,D] = eig(full(PA), full(PB));
        % Find the ones we're looking for.
        N = disc.dimension;
        idx = nearest(diag(D), P*V, sigma, min(k, N), N,disc);
        V = V(:, idx);
        D = D(idx, idx);
    else
        % FIXME: Experimental.
        [V, D] = eigs(PA, PB, k, sigma);
    end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Returns index vector that sorts eigenvalues by the given criterion.
function idx = nearest(lam, V, sigma, k, N, disc)

if ( isnumeric(sigma) )
    if ( isinf(sigma) )
        [junk, idx] = sort(abs(lam), 'descend');
    else
        [junk, idx] = sort(abs(lam-sigma));
    end
else
    switch upper(sigma)
        case 'LR'
            [junk, idx] = sort(real(lam), 'descend');
        case 'SR'
            [junk, idx] = sort(real(lam));
        case 'LI'
            [junk, idx] = sort(imag(lam), 'descend');
        case 'SI'
            [junk, idx] = sort(imag(lam));
        case 'LM'
            [junk, idx] = sort(abs(lam), 'descend');
            % case 'SM' already converted to sigma = 0
        otherwise
            error('CHEBFUN:linop:eigs:sigma', 'Unidentified input ''sigma''.');
    end
end

% Delete infinite values. These can arise from rank deficiencies in the
% RHS matrix of the generalized eigenproblem.
idx( ~isfinite(lam(idx)) ) = [];

% Propose to keep these modes.
queue = 1:min(k, length(idx));
keeper = false(size(idx));
keeper(queue) = true;

%%
% Screen out spurious modes. These are dominated by high frequency for all
% values of N. (Known to arise for some formulations in generalized
% eigenproblems, specifically Orr-Sommerfeld.)

% Grab some indices
tenPercent = ceil(N/10); % We are the 10%
iif10 = 1:tenPercent;    % Indices of first 10%
ii90 = tenPercent:N;     % Indices of last 90%
ii10 = (N-tenPercent):N; % Indices of last 10%

% Check for high frequency energy (indicative of spurious eigenvalues) in
% each of the remaining valid eigenfunctions.
isFun = disc.source.isFunVariable;
while ~isempty(queue)
    j = queue(1);

    vcoeff = mat2poly(disc, V(:,idx(j)));
    vcoeff = vcoeff(isFun);
    vcoeffsq = 0;
    for i = 1:numel(vcoeff)
        for q = 1:numel(vcoeff{i})
            newcoeff2 = vcoeff{i}{q}.*conj(vcoeff{i}{q});
            lnc2 = length(newcoeff2);
            lvcs = length(vcoeffsq);
            if ( lnc2 > lvcs )
                % Pad with leading zeros
                vcoeffsq = [ zeros(lnc2-lvcs,1); vcoeffsq ];
            else
                % Only the most significant rows affected
                rows = lvcs-lnc2+1:lvcs;
                vcoeffsq(rows) = vcoeffsq(rows) + newcoeff2;
            end
        end
    end
    vcoeff = sqrt( flipud( sum(vcoeffsq,2) ) );

    % Recipe: More than half of the energy in the last 90% of the Chebyshev
    % modes is in the highest 10% modes, and the energy of the last 90% is
    % not really small (1e-8) compared to the first 10% (i.e. is not noise).
    norm90 = norm(vcoeff(ii90)); % Norm of last 90%
    norm10 = norm(vcoeff(ii10)); % Norm of last 10%
    normFirst10 = norm(vcoeff(iif10)); % Norm of first 10%
    if ( norm10 > 0.5*norm90 && norm90 > 1e-8*normFirst10 )
        keeper(j) = false;
        if queue(end) < length(idx)
            m = queue(end)+1;
            keeper(m) = true;
            queue = [queue(:); m];
        end
    end
    queue(1) = [];

end

%%

% Return the keepers.
idx = idx( keeper );

end

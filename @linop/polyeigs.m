function varargout = polyeigs(varargin)
%POLYEIGS    Polynomial eigenvalues and eigenfunctions of a linear operator.
%   Important (1): While you can construct a LINOP and apply this method, the
%   recommended procedure is to use CHEBOP/POLYEIGS instead.
%   Important (2): A CHEBOPPREF object PREFS has to be passed. When this method
%   is called via CHEBOP/POLYEIGS, PREFS is inherited from the CHEBOP level.
%
% [X,E] = POLYEIG(A0,A1,..,Ap,K) solves the polynomial eigenvalue problem
% of degree p:
%    (A0 + lambda*A1 + ... + lambda^p*Ap)*x = 0.
% The input is p+1 linops, A0, A1, ..., Ap and the output is an inf-by-K
% chebfun quasimatrix, X, whose columns are the K least oscillatory
% eigenfunctions, and a vector of length k, E, whose elements are the
% eigenvalues.
%    for j = 1:K
%       lambda = E(j)
%       u = X(:,j)
%       A0(u) + lambda*A1(u) + ... + lambda^p*Ap(u) %is approximately 0.
%    end
% Any boundary conditions, continuity constraints, or other constraints on
% the eigenvectors should be homogeneous and specified in A0. Inhomogenous
% constraints are set to zero, and constarints in A1, ..., Ap are ignored.
% K defaults to 6 if not specified.
%
% E = POLYEIGS(A0,A1,..,Ap,K) is a vector of length k whose elements are
% the K least oscillatory eigenvalues of the polynomial eigenvalue problem.
%
% EIGS(A0,A1,..,Ap,K,SIGMA) also finds K solutions to the polynomial
% eigenvalue problem. If SIGMA is a scalar, the eigenvalues found are the
% ones closest to SIGMA. Other possibilities are 'LR' and 'SR' for the
% eigenvalues of largest and smallest real part, and 'LM' (or Inf) and 'SM'
% for largest and smallest magnitude. SIGMA must be chosen appropriately
% for the given operator; for example, 'LM' for an unbounded operator will
% fail to converge!
%
% Similarly to LINOP/EIGS, this routine uses the built-in POLYEIG on dense
% matrices of increasing size, stopping when the targeted eigenfunctions
% appear to have converged, as determined by the chebfun constructor.
%
% Example:
%
% d = [-1 1];
% x = chebfun('x', d);
% A = linop( operatorBlock.diff(d, 2) );
% E = functionalBlock.eval(d);
% A = addbc(A, E(-1), 0);
% A = addbc(A, E(1), 0);
% B = -linop( operatorBlock.mult(x)*operatorBlock.diff(d) );
% C = linop( operatorBlock.eye() );
% prefs = cheboppref();
% prefs.discretization = @chebcolloc2;
% [V,D] = polyeigs(A, B, C, 6, prefs)

% Copyright 2022 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Parsing inputs.
A = {};
k = [];       % will be made default value below
sigma = [];   % default 'auto' mode
prefs = cheboppref();
gotk = false; % until we detect a value of k in inputs
for j = 1:nargin
    item = varargin{j};
    if ( isa(item, 'linop') )
        % Generalized operator term
        A{end+1} = item;
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

% Degree of the polynomial eigenvalue equation.
p = numel(A)-1;

% Check for unbounded domains:
if ( ~all(isfinite(A{1}.domain)) )
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
m = size(A{1}, 2);
if ( m ~= size(A{1}, 1) )
    error('CHEBFUN:LINOP:eigs:notSquare','Block size must be square.')
end

% The domains must be merged before deriving the continuity equations for A:
dom = cellfun(@(A) A.domain, A, 'uniformoutput', false);
dom = domain.merge(dom{:});
for j = 1:p+1
    A{j}.domain = dom;
end

% Set up the discretization of A:
if ( isa(discType, 'function_handle') )
    % Create a discretization object
    discA = discType(A{1});

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
constructor = str2func( class(discA) );   % constructor handle.
discA = {discA};
for j = 2:p+1
    discA{j} = constructor(A{j});
end

% We can ignore constraints and continuity--enforced on A[0].
for j = 2:p+1
    if ( ~isempty(discA{j}.source.constraint) )
        discA{j}.source.constraint = [];
        warning('CHEBFUN:LINOP:eigs:constraints', ...
                'Constraints on A[%d] are ignored.', j-1);
    end
    if ( ~isempty(discA{j}.source.continuity) )
        discA{j}.source.continuity = [];
        warning('CHEBFUN:LINOP:eigs:continuity', ...
                'Continuity conditions on A[%d] are ignored.', j-1)
    end       
end

% Merge the all discretizations:
[discA{1:p+1}] = merge(discA{:});

if ( isempty(A{1}.continuity) )
     % Apply continuity conditions:
     discA{1}.source = deriveContinuity(discA{1}.source);
end

% 'SM' is equivalent to eigenvalues nearest zero.
if ( strcmpi(sigma, 'SM') )
    sigma = 0;
end

% Information required for finding the eigenvalues and functions.
numInts = discA{1}.numIntervals;
isFun = isFunVariable(A{1});

% Automatic mode: find the best sigma by going where the convergence appears to
% be fastest.
if ( isempty(sigma) )
    % Try to determine where the 'most interesting' eigenvalue is.
    discA{1}.dimension = 33*ones(1, numInts);
    [V1, D1] = getEigenvalues(discA, 33, 0);
    discA{1}.dimension(:) = 65;
    [V2, D2, P] = getEigenvalues(discA, 33, 0);
    lam1 = diag(D1);
    lam2 = diag(D2);
    dif = bsxfun(@minus, lam1.', lam2);
    delta = min( abs(dif) );   % diffs from 33->65
    % These are the significantly large differences from 33->65.
    bigDel = (delta > 1e-12*norm(lam1,Inf));

    % Trim off things that are still changing a lot (relative to new size).
    lam1b = lam1;
    lam1b(bigDel) = 0;
    bigDel = logical((delta > 1e-3*norm(lam1b, inf)) + bigDel);

    if ( all(bigDel) )
        % All values changed somewhat-- choose the one changing the least.
        [ignored, idx] = min(delta);
        sigma = lam1(idx);
    else
        % One by one, convert the eigenvectors to functions and check their cheb
        % expansion coefficients.
        U = partition(discA{1}, P*V2);  % each cell is array valued, for one variable

        % Combine the different variable components into a single variable for
        % coefficient conversion.
        Z = 0;
        for j = ( find(isFun) )
            Z = Z + U{j};
        end

        % Convert the discrete Z values to CHEBFUN
        z = toFunctionOut(discA{1}, Z);

        % Obtain all coefficients to use below
        coeffs = get(z, 'coeffs', 1);
        
        % Compute the 1-norm of the polynomial expansions, summing over smooth
        % pieces, for all columns.
        onenorm = 0;
        for j = 1:discA{1}.numIntervals
            onenorm = onenorm + sum(abs(coeffs{j}), 1 ).';
        end
        
        [ignored, index] = min(onenorm);
        sigma = lam2(index);
    end
end

% Linear combination coefficients for convergence test. The convergence of the
% combination is the same as the worst constituent function. The nontrivial
% coefficents are to make accidental cancellations extremely unlikely.
coeff = 1./(2*(1:k)');
for dim = [dimVals NaN]
    [V, D, P] = getEigenvalues(discA, k, sigma);

    % Combine the eigenfunctions into a composite.
    v = V*coeff(1:size(V,2));

    % Convert the different components into cells
    u = partition(discA{1}, P*v);

    % Test the happiness of the function pieces:
    vscale = zeros(1, sum(isFun));   % intrinsic scaling only
    [isDone, cutoff] = testConvergence(discA{1}, u(isFun), vscale, prefs);

    if ( all(isDone) )
        break
    elseif ( ~isnan(dim) )
        % Update the discretiztion dimension on unhappy pieces:
        discA{1}.dimension(~isDone) = dim;
    end

end

if ( ~isDone )
    warning('LINOP:EIGS:convergence', ...
        ['Maximimum dimension reached. Solution may not have converged.\n' ...
        'Please see help cheboppref.maxDimension for more details.']);
end

% Detect finite rank operators.
if ( size(D,1) < k )
    if ( gotk )
        warning('CHEBFUN:LINOP:eigs:rank',...
            'Input has finite rank, only %d eigenvalues returned.', size(D,1));
    end
    k = size(D,1);
end

% Sort eigenvalues:
[D, idx] = sort(D);
V = V(:,idx);

if ( nargout < 2 )  % Return the eigenvalues only
    varargout = { D };
else            % Unwrap the eigenvectors for output

    u = mat2fun(discA{1}, P*V, cutoff);

    % For normalizing eigenfunctions, so that they always have the same sign:
    signMat = [];
    
    % Find the norm in each eigenfunction (aggregated over variables).
    nrmsq = zeros(1,k);
    for j = 1:length(u)
        if ( isFun(j) )
            % Compress the representation.
            u{j} = simplify(u{j});
            if ( isempty(signMat) )
                % Find what domain we are working on:
                dom = domain(u{j});
                % Arbitrary point just to the right of the middle of the domain:
                fevalPoint = dom(1) + diff([dom(1) dom(end)])*.500023981;
                % Find out what sign the real part of the function have there:
                fevalSigns = sign(real(feval(u{j}, fevalPoint)));
                % Diagonal matrix with elements equal to the sign at our
                % arbitrary point. Add 0.1 and take signs again to ensure we
                % don't end up with any zeros (in case we were very unlucky).
                signMat = diag(sign(fevalSigns + 0.1));
            end
        end
        nrmsq = nrmsq + sum(u{j}.*conj(u{j}), 1);
    end
    
    % Normalize each eigenfunction.
    scale = diag( 1./sqrt(nrmsq') );
    for j = 1:length(u)
        u{j} = u{j}*scale*signMat;
    end

    % TODO: Can we move this to the CHEBMATRIX constructor?
    % NOTE: The following is required because block entries of a CHEBMATRIX
    % should only contain scalar objects (in particular, _not_ array-valued
    % CHEBFUNS or quasimatrices). Here we unwrap everything so that each
    % component of each eigenfunction is a single entry in a cell array.
    for j = 1:numel(u)
        % Convert each solution to it's own entry in a cell.
        u{j} = num2cell(u{j});
    end
    u = chebmatrix(vertcat(u{:}));

    % Output:
    varargout = {u, D};
end

end
% END OF MAIN FUNCTION

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [V, D, P] = getEigenvalues(discA, k, sigma)
% Formulate the discrete problem and solve for the eigenvalues

    p = numel(discA)-1;
    
    % Discretize the LHS operator (incl. constraints/continuity):
    [PA, P, C, ignored, PS] = matrix(discA{1});
    PA = {PA};
    for j = 2:p+1
        discA{j}.dimension = discA{1}.dimension;
        PB = matrix(discA{j});
        PA{j} = [ zeros(size(C)) ; PB ];
    end
    % Enforce dense matrices:
    PA = cellfun(@full, PA, 'uniformoutput', false);

    % Compute eigenvalues.
    [V, D] = polyeig(PA{:});
    lam = D;
    
    % Find the ones we're looking for.
    lam = deflate(lam, size(C,1));
    idx = nearest(lam, sigma);
    idx = filter(idx, P*V, k, discA{1});
    
    % Extract them:
    V = V(:,idx);
    D = D(idx);
        
end

function lam = deflate(lam, m)
% DEFLATE(LAM, M) forces that the M largest eigenvalues are deflated to INF.

[junk, idx] = sort(abs(lam), 'descend');
lam(idx(1:m)) = inf;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function idx = nearest(lam, sigma)
% Returns index vector that sorts eigenvalues by the given criterion.

if ( isnumeric(sigma) )
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

function idx = filter(idx, V, k, disc)
% Screen out spurious modes. These are dominated by high frequency for all
% values of N. (Known to arise for some formulations in generalized
% eigenproblems, specifically Orr-Sommerfeld.)

% TODO: Explain this in more detail.

N = disc.dimension;
k = min(k, N);

% Propose to keep these modes.
queue = 1:min(k, length(idx));
keeper = false(size(idx));
keeper(queue) = true;

% Grab some indices
tenPercent = ceil(N/10); % We are the 10%
iif10 = 1:tenPercent;    % Indices of first 10%
ii90 = tenPercent:N;     % Indices of last 90%
ii10 = (N-tenPercent):N; % Indices of last 10%

% Check for high frequency energy (indicative of spurious eigenvalues) in
% each of the remaining valid eigenfunctions.
isFun = disc.source.isFunVariable;
while ( ~isempty(queue) )
    j = queue(1);

    vcoeff = mat2poly(disc, V(:,idx(j)));
    vcoeff = vcoeff(isFun);
    vcoeffsq = 0;
    for i = 1:numel(vcoeff)
        for q = 1:numel(vcoeff{i})
            % TODO: The flipud below is required to make sure that the 
            % algorithm, designed for the old ordering of cheb-coeffs, continues
            % to work. One can remove the following flipud but then carefull 
            % changes will be needed in this function.
            vcoeff{i}{q} = flipud(vcoeff{i}{q});
            newcoeff2 = vcoeff{i}{q}.*conj(vcoeff{i}{q});
            lnc2 = length(newcoeff2);
            lvcs = length(vcoeffsq);
            if ( lnc2 > lvcs )
                % Pad with leading zeros
                vcoeffsq = [ zeros(lnc2 - lvcs,1) ; vcoeffsq ; ]; %#ok<AGROW>
                lvcs = length(vcoeffsq);
            end
            % Only the most signifi        256         512cant rows affected
            rows = (lvcs - lnc2 + 1):lvcs;
            vcoeffsq(rows) = vcoeffsq(rows) + newcoeff2; %#ok<AGROW>
        end
    end
    vcoeff = sqrt( flipud(sum(vcoeffsq, 2)) );

    % Recipe: More than half of the energy in the last 90% of the Chebyshev
    % modes is in the highest 10% modes, and the energy of the last 90% is
    % not really small (1e-8) compared to the first 10% (i.e. is not noise).
    norm90 = norm(vcoeff(ii90)); % Norm of last 90%
    norm10 = norm(vcoeff(ii10)); % Norm of last 10%
    normFirst10 = norm(vcoeff(iif10)); % Norm of first 10%
    if ( norm10 > 0.5*norm90 && norm90 > 1e-8*normFirst10 )
        keeper(j) = false;
        if queue(end) < length(idx)
            m = queue(end) + 1;
            keeper(m) = true;
            queue = [queue(:); m];
        end
    end
    queue(1) = [];
    
end

% Return the keepers.
idx = idx( keeper );

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

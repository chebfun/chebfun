function [v, S] = null(A, prefs, nullity)
%NULL   Null space of a LINOP.
%   Important (1): While you can construct a LINOP and apply this method, the
%   recommended procedure is to use CHEBOP/NULL instead.
%   Important (2): A CHEBOPPREF object PREFS has to be passed. When this method
%   is called via CHEBOP/NULL, PREFS is inherited from the CHEBOP level.
%
%   Z = NULL(A, PREFS) returns a CHEBMATRIX with orthonormal columns which span 
%   the null space of the LINOP A. That is, A*Z has negligible elements, 
%   SIZE(Z, 2) is the nullity of A, and Z'*Z = I. A may contain linear 
%   boundary conditions, but they will be treated as homogeneous. The nullity 
%   is determined by comparing the differential order of the system and the
%   number of supplied boundary conditions.
%
%   Z = NULL(A, PREFS, K) or Z = NULL(A, K) attempts to find K null vectors. If
%   the nullity of A is determined to be less than K then an warning is thrown.
%   This is useful in situations where the nullity is known in advance and the
%   algorithm struggles to determine it automatically.
%
% Example:
%   d = domain(0, pi);
%   A = diff(d);
%   prefs = cheboppref();
%   prefs.discretization = @chebcolloc2;
%   V = null(A, prefs);
%   norm(A*V)
%
% See also LINOP/SVDS, LINOP/EIGS, NULL.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Discretization type:
discType = prefs.discretization;
% Take the standard BVP tolerance:
tol = prefs.bvpTol;

if ( nargin < 3 )
    nullity = [];
elseif ( nullity == 0 )
    nullity = [];
end

% Check for square operator. (This is not strict enough, technically.)
m = size(A, 2);
if ( m ~= size(A, 1) )
    error('CHEBFUN:LINOP:eigs:notSquare', 'Block size must be square.')
end

% Set up the discretization:
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

if ( isempty(A.continuity) )
     % Apply continuity conditions:
     discA.source = deriveContinuity(discA.source);
end

if ( isempty(nullity) )
    nullity = sum(discA.projOrder) - ...
    ( size(discA.source.continuity.values, 1) + ...
      size(discA.source.constraint.values, 1) );
end

% Boundary conditions are not applied, so we want square operators:
discA.dimAdjust = 0*discA.dimAdjust;
discA.projOrder = 0*discA.projOrder;

% Information required for finding the nullity and null space.
isFun = isFunVariable(A);

% Linear combination coefficients for convergence test. The convergence of the
% combination is the same as the worst constituent function. The nontrivial
% coefficents are to make accidental cancellations extremely unlikely.
coeff = @(n) 1./(2*(1:n).');

for dim = dimVals

    % Get discrete null vectors:
    [V, P, S] = getNullVectors(discA, tol, nullity);

    % Combine the singular vectors into a composite.
    v = V*coeff(size(V, 2));
    % Convert the different components into cells
    v = partition(discA, P*v);
    
    % Test the happiness of the function pieces:
    vscale = zeros(1, sum(isFun));  % intrinsic scaling only.
    isDone = testConvergence(discA, v(isFun), vscale, prefs);
    
    if ( all(isDone) )
        break
    else
        % Update the discretiztion dimension on unhappy pieces:
        discA.dimension(~isDone) = dim;
    end

end

if ( ~isempty(nullity) && max(S) > tol ) 
    warning('CHEBFUN:linop:null:nullity', ...
        ['Number of requested null-vectors exceeds computed ', ...
         'nullity of the operator.\n', ...
        '(Relative) singular value of ' num2str(norm(S, inf)) ' encountered.'])
    V = V(:,1:(sum(S<tol) - 1));
end

if ( isempty(V) )
    v = [];
    return
end

% Convert discrete data to functions:
v = mat2fun(discA, P*V);

% Simplify and orthogonalize:
if ( m == 1 )
    v{1} = qr(v{1});
    v{1} = simplify(v{1});
else % system of eqns
    [~, R] = qr(join(v{:}));
    for j = 1:numel(v)
        v{j} = mrdivide(v{j}, R, 'least-squares');
        v{j} = simplify(v{j});
    end
end

% TODO: Can we move this to the CHEBMATRIX constructor? NOTE: The following is
% required because block entries of a CHEBMATRIX should only contain scalar
% objects (in particular, _not_ array-valued CHEBFUNS or quasimatrices). Here we
% unwrap everything so that each component of each eigenfunction is a single
% entry in a cell array.
for j = 1:numel(v)
    % Convert each solution to it's own entry in a cell.
    v{j} = num2cell(v{j});
end
v = chebmatrix(vertcat(v{:}));

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [V, P, S] = getNullVectors(discA, tol, nullity)
% Formulate the discrete problem and solve for the eigenvalues

    % Discretize the LHS operator (incl. constraints/continuity):
    [~, P, B, A] = matrix(discA);
    
    % Construct one big matrix from the unprojected block entries:
    A = cell2mat(A);
    
    % Compute the discrete SVD. (Note: It saves no time to call the built-in
    % NULL() method, and this simply calls the built-in SVD method.)
    [~, S, V] = svd(full(A), 0);
    S = diag(S)/S(1);

    % Numerical nullity:
    if ( isempty(nullity) )
        idx = sum(S < tol);
        S = S(end-idx+1:end);
    else
        idx = nullity + size(B, 1);
        S = S(end-idx+1:end);
    end

    % Extract null vectors:
    if ( idx ~= 0 )
        V = V(:,end+1-idx:end);        % Numerical null vectors.
        % Enforce additional boundary conditions:
        if ( ~isempty(B) )
            V = V*null(B*V);           % Store output in V.
        end
    else
        V = V(:,end); % Check for convergence in smallest singular value.
    end

end

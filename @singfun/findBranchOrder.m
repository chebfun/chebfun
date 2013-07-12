function branchOrder = findBranchOrder( op, singFlag )
% One test for blowup - perhaps less robust than test B,
% but more accurate when it works.  This is based on
% a test of monotonicity of the function and its derivative.
% We return 1 if there's a blowup, 0 if not.

% First get an estimate of the exponent
% by the pole order finder. These will be
% passed on to the branchOrderFinder as 
% upperbounds for the singularity order.
poleBound = findPoleOrder( op, singFlag );

branchOrder = zeros(1,2);
% distance of sample points from the end points
x = eps*(11:-1:2)';
% if a pole is expected at x = 1
if ( isnan(singFlag(2)) )
    fvalsRight = op(1-x);
    branchOrder(2) = branchOrderFinder( fvalsRight, x, poleBound(2) );
else
    % no singularity
    branchOrder(2) = 0;
end
% if a pole is expected at x = -1
if ( isnan(singFlag(1)) )
    fvalsLeft = op(-1+x);
    branchOrder(1) = branchOrderFinder( fvalsLeft, x, poleBound(1) );
else
    % no singularity
    branchOrder(1) = 0;
end
end

function branchOrder = branchOrderFinder( fvals, x, poleBound )
% Iteratively increase the proposed value. If
% the singularity at the right endpoint can be mollified, we can declare
% success.

branchOrder = poleBound-1;
% decimal search
tol=1e-12; % what should this be [TODO]
maxIter = 100;
% initial grid of size n: assume the order of the 
% singularity between poleOrder-1 and poleOrder
n = 10;
exponentGrid = linspace(branchOrder, branchOrder+1, n);
absFvals = abs(fvals);
smoothVals = absFvals;
nIter = 0;
while( abs(exponentGrid(end)-exponentGrid(1)) > tol && nIter <= maxIter )
    k = 0;
    while( all(diff(diff(smoothVals)) > 0) && k < n )
        k = k + 1;
        smoothVals = absFvals.*x.^exponentGrid(k);
    end
    if( k == n && all(diff(diff(smoothVals)) > 0) )
        % tried all exponents but failed
        warning( 'something' );
        branchOrder = poleBound;
    else
        % succeeded for some k, update the estimate and refine the grid
        smoothVals = absFvals;
        branchOrder = exponentGrid(k);
        exponentGrid = linspace( exponentGrid(k-1), exponentGrid(k), n );
        nIter = nIter + 1;
    end
end
if ( nIter >= maxIter )
    warning( 'da da da' );
    branchOrder = poleBound;
end
end
function singOrder = findSingOrder( op, singEnd )
% One test for blowup - perhaps less robust than test B,
% but more accurate when it works.  This is based on
% a test of monotonicity of the function and its derivative.
% We return 1 if there's a blowup, 0 if not.

% First get an estimate of the exponent
% by the pole order finder. These will be
% passed on to the branchOrderFinder as 
% upperbounds for the singularity order.
poleBound = singfun.findPoleOrder( op, singEnd );

% distance of sample points from the end points
x = eps*(11:-1:2)';
% if a pole is expected at x = 1
if ( strcmpi(singEnd, 'right') )
    fvalsRight = op(1-x);
    singOrder = singOrderFinder( fvalsRight, x, poleBound);
else if ( strcmpi(singEnd, 'left') )
        % if a pole is expected at x = -1
        fvalsLeft = op(-1+x);
        singOrder = singOrderFinder( fvalsLeft, x, poleBound);
    else
        error('CHEBFUN:SINGFUN:findSingOrder:unknownPref',...
                    'Blowup preference "%s" unknown', singEnd )
    end
end

end

function singOrder = singOrderFinder( fvals, x, poleBound )
% Iteratively increase the proposed value. If
% the singularity at the right endpoint can be mollified, we can declare
% success.

singOrder = poleBound-1;
% decimal search
maxIter = 100;
% initial grid of size n: assume the order of the 
% singularity between poleOrder-1 and poleOrder
n = 10;
exponentGrid = linspace(singOrder, singOrder+1, n);
absFvals = abs(fvals);
smoothVals = absFvals;
nIter = 0;
tol = singfun.pref.singfun.eps;
while( abs(exponentGrid(end)-exponentGrid(1)) > 100*tol && nIter <= maxIter )
    k = 0;
    while( all(diff(diff(smoothVals)) > 0) && k < n )
        k = k + 1;
        smoothVals = absFvals.*x.^exponentGrid(k);
    end
    if( k == n && all(diff(diff(smoothVals)) > 0) )
        % tried all exponents but failed
        singOrder = poleBound;
        return;
    else if( k == 0 )
            singOrder = poleBound;
            return;
        else
            % succeeded for some k, update the estimate and refine the grid
            smoothVals = absFvals;
            singOrder = exponentGrid(k);
            exponentGrid = linspace( exponentGrid(k-1), exponentGrid(k), n );
            nIter = nIter + 1;
        end
    end
end
end
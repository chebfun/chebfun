function singOrder = findSingOrder(op, singEnd)
%FINDSINGORDER   Finds the order of the algebraic singularity (a fractional
%   pole) in the function handle OP at x = 1 or -1 depending upon the string
%   'left' or 'right' passed in SINGEND. The method also works for poles, 
%   i.e. if the order of the singularity is an integer.
%    Example:
%   p = singfun.findSingOrder(@(x) 1./(1-x).^1.5, 'right' )
%   p = singfun.findSingOrder(@(x) 1./(1+x).^2.5, 'left' )
%   p = singfun.findSingOrder(@(x) 1./(1+x).^5, 'left' )
%
% See also FINDPOLEORDER and FINDSINGEXPONENTS

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

%%
% First get an estimate of the exponent by the pole order finder. This 
% will be passed on as upperbound for the singularity order.
poleBound = -singfun.findPoleOrder(op, singEnd );

% distance of sample points from the end points
x = eps*(11:-1:2)';
% if a pole is expected at x = 1
if ( strcmpi(singEnd, 'right') )
    fvalsRight = op(1 - x);
    singOrder = singOrderFinder( fvalsRight, x, poleBound);
else if ( strcmpi(singEnd, 'left') )
        % if a pole is expected at x = -1
        fvalsLeft = op(-1 + x);
        singOrder = singOrderFinder( fvalsLeft, x, poleBound);
    else
        error('CHEBFUN:SINGFUN:findSingOrder:unknownPref', ...
              'Blowup preference "%s" unknown', singEnd )
    end
end
% singOrder is a negative number.
singOrder = -singOrder;

% The algorithm does not support positive exponents >= 1. Such results are 
% garbage, so discard them. This is consistent with Chebfun V4. From a
% practial point of view, this is not a problem since a Chebfun with a
% barnch singularity of order > 1 converges. For such functions, the 
% cheb-coefficient decay is algebraic but strong enough to give us a reliable 
% representation.
if ( singOrder >= 1 )
    singOrder = 0;
end

end

function singOrder = singOrderFinder( fvals, x, poleBound )
%SINGORDERFINDER   Finds the order of the singularity based on function 
%   values FVALS given at (1-X). POLEBOUND is an integer upperbound 
%   of the singularity order.

singOrder = poleBound - 1;

% Try fractional exponents in the interval [poleBound-1, poleBound]
% Take an initial grid of size n: 
n = 11;
exponentGrid = linspace(poleBound - 1, poleBound, n);

% some initialisations which will be used later
absFvals = abs(fvals);
nIter = 0;

%%
% Tolerance for exponents
tol = singfun.pref.singfun.exponentTol;

% maximum number of iterations allowed
maxIter = 100; 
% A factor by 10 refinement algorithm for zooming in on the required
% fractional singularity order.
while( (abs(exponentGrid(end) - exponentGrid(1)) > tol) && (nIter <= maxIter) )
    k = 1;
    smoothVals = absFvals.*x.^exponentGrid(k);

    % This test for blowup is perhaps less robust than the test given in
    % SINGFUN.FINDPOLEORDER() but more accurate when it works. This is based 
    % on determining the convexity of the function via second order
    % differences.
    while( all(diff(diff(smoothVals)) > 0) && (k < n) )
        k = k + 1;
        % Try the next fractional exponent
        smoothVals = absFvals.*x.^exponentGrid(k);
    end
    
    if( (k == n) && all(diff(diff(smoothVals)) > 0) )
        % tried all exponents but failed
        singOrder = exponentGrid(n);
        return;
    else
        % update the estimate for exponent
        singOrder = exponentGrid(k);        
        if( k == 1 )
            % no further refinement possible.            
            return;
        else           
            % refine the grid                        
            exponentGrid = linspace(exponentGrid(k - 1), exponentGrid(k), n);
            nIter = nIter + 1;
        end
    end
end

end

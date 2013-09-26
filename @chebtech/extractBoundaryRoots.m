function [f, rootsLeft, rootsRight] = extractBoundaryRoots(f)

% Copyright 2013 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% Tolerance for a root:
tol = 500*eps*get(f, 'vscale');

% Values at ends:
endValues = abs([get(f, 'lval'), get(f, 'rval')]);

% Initialise the roots:
rootsLeft = 0;
rootsRight = 0;

% If there are no roots, there is nothing to do!
if ( all(endValues > tol) )
    return
end

c = f.coeffs;

while ( any(endValues < tol) )
    
    if ( endValues(1) < tol )
        % Root at the left.
        sgn = 1;   
        rootsLeft = rootsLeft + 1;
    else
        % Root at the right.
        sgn = -1;
        rootsRight = rootsRight + 1;
    end
    
    % Construct the matrix for the recurrence
    n = length(f);
    e = ones(n-1, 1);
    D = spdiags([.5*e, sgn*e, .5*e], 0:2, n-1, n-1); 
    D(1) = 1; %#ok<SPRIX>
    
    % Compute the new coefficients
    c = sgn*flipud(D\c(end-1:-1:1,:));
   
    % Construct new f
    f.values = f.chebpolyval(c);
    f.coeffs = c;
    
    % Update endValues:
    endValues = abs([get(f, 'lval'), get(f, 'rval')]);

end
   
end
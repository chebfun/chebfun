function [f, rootsLeft, rootsRight] = extractBoundaryRoots(f)
%EXTRACTBOUNDARYROOTS   Extract roots at the boundary points -1 and 1.
%   [F, ROOTSLEFT, ROOTSRIGHT] = EXTRACTBOUNDARYROOTS(F) returns a CHEBTECH G
%   which is free of roots at the boundary points -1 and 1. The multiplicity of
%   the boundary roots at -1 and 1 are ROOTSLEFT and ROOTRIGHT respectively.
%
% See also ROOTS.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Grab the size of F:
m = size(f, 2);

% Tolerance for a root (we will loosen this with each run of the loop below if
% there are multiple roots):
tol = 1e2*f.vscale*f.epslevel;

% Values at ends:
endValues = abs([feval(f, -1); feval(f, 1)]);

% Initialise the multiplicity of the roots:
rootsLeft = zeros(1, m);
rootsRight = zeros(1, m);

% If there are no roots, there is nothing to do!
if ( all(min(endValues, [], 1) > tol) )
    return
end

% Grab the coefficients of F:
c = f.coeffs;

while ( any(min(endValues, [], 1) <= tol) )
    
    if ( any(endValues(1, :) <= tol) )
        % Root at the left.
        sgn = 1;
        ind = find(endValues(1, :) <= tol);
        rootsLeft(ind) = rootsLeft(ind) + 1;
    else
        % Root at the right.
        sgn = -1;
        ind = find(endValues(2, :) <= tol);
        rootsRight(ind) = rootsRight(ind) + 1;
    end
    
    % Construct the matrix for the recurrence:
    n = length(f);
    e = ones(n-1, 1);
    D = spdiags([.5*e, sgn*e, .5*e], 0:2, n-1, n-1);
    D(1) = 1; %#ok<SPRIX>
    
    % Compute the new coefficients:
    c(2:end,ind) = sgn*flipud(D\c(end-1:-1:1,ind));
    
    % Pad zero at the highest coefficients:
    c(1,ind) = 0; 
    
    % Update the coefficients.  Note that we don't need to update the values 
    % here, since only coefficients are used in this while loop.  (feval(), 
    % below, calls CLENSHAW, which only uses coefficients.)
    f.coeffs = c;
    
    % Update endValues:
    endValues = abs([feval(f, -1); feval(f, 1)]);
    
    % Loosen the tolerance for checking multiple roots:
    tol = 1e2*tol;
    
end

% Call simplify to simplify and update the vscale:
f = simplify(f);

end
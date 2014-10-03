function [f, rootsLeft, rootsRight] = extractBoundaryRoots(f, numRoots)
%EXTRACTBOUNDARYROOTS   Extract roots at the boundary points -1 and 1.
%   [G, ROOTSLEFT, ROOTSRIGHT] = EXTRACTBOUNDARYROOTS(F) returns a CHEBTECH G
%   which is free of roots at the boundary points -1 and 1. The multiplicity of
%   the boundary roots at -1 and 1 are ROOTSLEFT and ROOTRIGHT respectively.
%
%   [G, ROOTSLEFT, ROOTSRIGHT] = EXTRACTBOUNDARYROOTS(F, NUMROOTS) returns a 
%   CHEBTECH G which has been peeled off roots from the boundaries with the
%   multiplicities specified by ROOTSLEFT and ROOTSRIGHT for the left boundary 
%   and the right boundary, respectively. NUMROOTS is a 2xM matrix with M being 
%   the number of the columns of the array-valued CHEBTECHs specifying the 
%   multiplicities of the boundary roots that EXTRACTBOUNDARYROOTS aims to 
%   extract. The first row and the second row of NUMROOTS correspond to the left
%   and the right boundary, respectively. Note that there is no warning thrown,
%   when ROOTSLEFT and ROOTSRIGHT are not consistent with NUMROOTS.
%
% See also ROOTS.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( nargin == 1 )
    numRoots = [];
end

% Grab the size of F:
m = size(f, 2);

% Initialise the multiplicity of the roots:
rootsLeft = zeros(1, m);
rootsRight = zeros(1, m);

% Tolerance for a root (we will loosen this with each run of the loop below if
% there are multiple roots):
tol = 1e3*f.vscale.*f.epslevel;

% Values at ends:
endValues = abs([feval(f, -1); feval(f, 1)]);

% If there are no roots, there is nothing to do!
if ( all(min(endValues, [], 1) > tol) )
    return
end

% Grab the coefficients of F:
c = f.coeffs;

while ( ( ( nargin == 1 ) && any( min(endValues, [], 1) <= tol ) ) ...
        || ( ( nargin > 1 ) && any( any( numRoots > 0 ) ) ) )
    
    if ( nargin == 1 )
        if ( any(endValues(1, :) <= tol) )
            % Root at the left:
            sgn = 1;
            ind = find(endValues(1, :) <= tol);
            rootsLeft(ind) = rootsLeft(ind) + 1;
        else
            % Root at the right:
            sgn = -1;
            ind = find(endValues(2, :) <= tol);
            rootsRight(ind) = rootsRight(ind) + 1;
        end
    else
        if ( any( numRoots(1,:) ) )
            % Root at the left:
            ind  = ( endValues(1, :) <= tol );
            indNumRoots = ( numRoots(1,:) > 0 );
            if ( ind == indNumRoots )
                sgn = 1;
                numRoots(1, ind) = numRoots(1, ind) - 1;
                rootsLeft = rootsLeft + 1;
            else
                numRoots(1,:) = 0;
                continue
            end
            
        elseif ( any( numRoots(2,:) ) )
            % Root at the right:
            ind  = ( endValues(2, :) <= tol );
            indNumRoots = ( numRoots(2,:) > 0 );
            if ( ind == indNumRoots )
                sgn = -1;
                numRoots(2, ind) = numRoots(2, ind) - 1;
                rootsRight = rootsRight + 1;
            else
                numRoots(2,:) = 0;
                continue
            end
        end
    end
    
    % Construct the matrix for the recurrence:
    n = length(f);
    e = ones(n-1, 1);
    D = spdiags([.5*e, sgn*e, .5*e], 0:2, n-1, n-1);
    D(1) = 1; %#ok<SPRIX>
    
    % Compute the new coefficients:
    c(1:end-1,ind) = sgn*(D\c(2:end,ind));
    
    % Pad zero at the highest coefficients:
    c(end,ind) = 0; 
    
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

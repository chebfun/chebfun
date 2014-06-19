function out = isdecay(f)
%ISDECAY   Test if a CHEBTECH decays faster than a single root at endpoints.
%   ISDECAY(F) returns a 1x2 row vector each of which indicates whether F
%   vanishes at one of the endpoints faster than a single root. An entry TRUE is
%   returned if F has a boundary root with multiplicity larger than one, FALSE
%   otherwise. 
%
%   Note that ISDECAY is designed for and expected to be called only by UNBNDFUN
%   class for handling functions defined on unbounded domains.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

out = zeros(2, size(f, 2));

% Set a tolerance:
tol = 1e2*get(f, 'epslevel').*get(f, 'vscale');

%% If F is a constant:
if ( length(f) == 1 )
    mask = ( f.coeffs < tol ) | ( f.coeffs == 0 );
    if ( any(mask) )
        out(:, mask) = 1;
    end
    return
end

%% Left endpoint:
rootsLeft = zeros(2, size(f, 2));
maskLeft = abs(get(f, 'lval')) < tol;
if ( any(maskLeft) )
    rootsLeft(1, maskLeft) = 1;
    
    % Extract a single boundary root at the left end points:
    g = extractBoundaryRoots(f, rootsLeft);
    
    % Check decaying speed:
    mask = abs(get(g, 'lval')) < 1e3*tol;
    if ( any(mask) )
        out(1, mask) = 1;
    end
end

%% Right endpoint:
rootsRight = zeros(2, size(f, 2));
maskRight = abs(get(f, 'rval')) < tol;
if ( any(maskRight) )
    rootsRight(2, maskRight) = 1;
    
    % Extract a single boundary root at the left end points:
    g = extractBoundaryRoots(f, rootsRight);
    
    % Check decaying speed:
    mask = abs(get(g, 'rval')) < 1e3*tol;
    if ( any(mask) )
        out(2, mask) = 1;
    end
end

end
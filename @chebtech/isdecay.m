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

out = zeros(1, 2);

% Set a tolerance:
tol = 1e2*get(f, 'epslevel')*get(f, 'vscale');

%% If G is a constant:

if ( length(f) == 1 )
    if ( ( f.coeffs < tol ) || ( f.coeffs == 0 ) )
        out = ones(1, 2);
    end
    return
end

%% Left endpoint:

if ( abs(get(f, 'lval')) < tol )
    % Extract a single boundary root at the left end points:
    g = extractBoundaryRoots(f, [1 ; 0]);
    
    % Check decaying speed:
    if ( abs(get(g, 'lval')) < 1e3*tol )
        out(1) = 1;
    end
end

%% Right endpoint:

if ( abs(get(f, 'rval')) < tol )
    % Extract a single boundary root at the left end points:
    g = extractBoundaryRoots(f, [0 ; 1]);
    
    % Check decaying speed:
    if ( abs(get(g, 'rval')) < 1e3*tol )
        out(2) = 1;
    end
end

end

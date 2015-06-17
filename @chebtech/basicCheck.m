function [ishappy, epslevel, cutOff] = basicCheck(f, varargin)
%BASICCHECK   Attempt to trim trailing Chebyshev coefficients in a CHEBTECH.
%   [ISHAPPY, EPSLEVEL, CUTOFF] = BASICCHECK(F, VALUES, VSCL, PREF) calls 
%   PLATEAUCHECK to check for convergence and always returns 
%   EPSLEVEL = eps and CUTOFF = LENGTH(F).
%
% See also LINOPV4CHECK, STRICTCHECK, CLASSICCHECK.

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Check convergence using plateauCheck:
[ishappy, epslevel, cutOff] = plateauCheck(f, varargin{:});

% do not trim coefficients
cutOff = length(f);

% set epslevel
epslevel = eps + 0*f.vscale;

end

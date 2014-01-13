function f = power(f, b)
% .^   SINGFUN power.
%   F.^G returns a SINGFUN F to the scalar power G, a scalar F to the SINGFUN
%   power G, or a SINGFUN F to the SINGFUN power G. F and or G may be complex. 
%   Note that it is assumed that F is non-zero on its domain. If F has zeros, 
%   then the output is garbage without throwing a warning.
%
%   H = POWER(F, G) is called for the syntax 'F .^ G'.
%
% See also SQRT.

% Copyright 2013 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% Extract roots from the boundaries and incrememnt the exponents accordingly:
f = extractBoundaryRoots(f);

% Modify the exponents:
f.exponents = b*f.exponents;

% Call POWER@SMOOTHFUN to update f.smoothPart (the output of which is expected 
% to be smooth):
f.smoothPart = power(f.smoothPart, b);

end
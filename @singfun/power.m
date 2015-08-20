function f = power(f, b)
%.^   SINGFUN power.
%   F.^G returns a SINGFUN F to the scalar power G, a scalar F to the SINGFUN
%   power G, or a SINGFUN F to the SINGFUN power G. F and or G may be complex. 
%
%   This function assumes that the curve traced out by F in the complex plane
%   both (1) does not come too close to zero except at the domain boundaries 
%   +/- 1 and (2) does not cross over the branch cut in POWER along the negative
%   real axis.  That is, F should not vanish at any point of (-1, 1), and the
%   imaginary part of F should not vanish at any point of (-1, 1) where the real
%   part of F is negative.  If any of these assumptions are violated, garbage
%   may be returned with no warning.
%
%   H = POWER(F, G) is called for the syntax 'F .^ G'.
%
% See also SQRT.

% Copyright 2015 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% Extract roots from the boundaries and incrememnt the exponents accordingly:
f = extractBoundaryRoots(f);

% Modify the exponents:
f.exponents = b*f.exponents;

% Call POWER@SMOOTHFUN to update f.smoothPart (the output of which is expected 
% to be smooth):
f.smoothPart = power(f.smoothPart, b);

% Simplify the exponents:
f = simplifyExponents(f);

end

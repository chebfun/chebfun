function out = isequal(f, g)
%ISEQUAL   Test if SINGFUN objects F and G are equal.
%   ISEQUAL(F, G) returns TRUE if the SINGFUN objects F and G have the same
%   underlying SMOOTHPART and the same EXPONENTS. By same EXPONENTS we mean
%   that they agree up to SINGFUN tolerance.

% Copyright 2015 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

tol = chebfunpref().blowupPrefs.exponentTol;

if ( isa(f, 'smoothfun') )
    f = singfun.smoothFun2SingFun(f);
end

if ( isa(g, 'smoothfun') )
    g = singfun.smoothFun2SingFun(g);
end

out = all(abs(f.exponents - g.exponents) < tol) && ...
    isequal(f.smoothPart, g.smoothPart);

end

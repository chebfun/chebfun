function f = flipud(f)
%FLIPUD   Flip/reverse a BNDFUN.
%   G = FLIPUD(F), where the BNDFUN F has a domain [a,b], returns a BNDFUN G
%   defined on [a,b], such that G(x) = F(-x) for all x in [a,b].

% Copyright 2013 by The University of Oxford and The Chebfun Developers. 
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.

% Flip the onefun:
f.onefun = flipud(f.onefun);

end
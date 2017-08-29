function f = flipud(f)
%FLIPUD   Flip/reverse a BNDFUN.
%   G = FLIPUD(F), where the BNDFUN F has a domain [a,b], returns a BNDFUN G
%   defined on [a,b], such that G(x) = F(a+b-x) for all x in [a,b], i.e., FLIPUD
%   flips the BNDFUN around the center of the domain it is defined on.

% Copyright 2017 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% Flip the onefun:
f.onefun = flipud(f.onefun);

end

function H = xor(F,G)
%XOR   Logical CHEBFUN EXCLUSIVE OR.
%     XOR(S,T) is the logical symmetric difference of CHEBFUNs S and T.
%     The result is logical 1 (TRUE) where either S or T, but not both, is
%     nonzero.  The result is logical 0 (FALSE) where S and T are both zero
%     or nonzero.  S and T must have the same dimensions (or one can be a
%     scalar).

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

H = or(F, G) & ~and(F, G);

end

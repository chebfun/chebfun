function H = xor(F,G)
% XOR Logical chebfun EXCLUSIVE OR.
%     XOR(S,T) is the logical symmetric difference of chebfuns S and T.
%     The result is logical 1 (TRUE) where either S or T, but not both, is
%     nonzero.  The result is logical 0 (FALSE) where S and T are both zero
%     or nonzero.  S and T must have the same dimensions (or one can be a
%     scalar).

H = or(F,G) & ~and(F,G);
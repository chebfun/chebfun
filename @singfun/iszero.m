function out = iszero(f)
%ISZERO   True for zero SINGFUN objects.
%   ISZERO(F) returns logical TRUE if the smooth part of  F is zero and FALSE
%   otherwise.

out = iszero(f.smoothPart);

end

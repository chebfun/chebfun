function out = iszero(f)
%ISZERO    True for zero FUN objects.
%   ISZERO(F) returns logical TRUE is F.onefun is identically zero and logical
%   FALSE otherwise.

% Call ISZERO() of the ONEFUN:
out = iszero(f.onefun);

end
function out = iszero(f)
%ISZERO    True for zero BNDFUN objects.
%   ISZERO(F) returns logical TRUE is F.onefun.values has only zero entries and 
%   logical FALSE otherwise.

out = iszero(f.onefun);

end
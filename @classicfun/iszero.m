function out = iszero(f)
%ISZERO    True for zero CLASSICFUN objects.
%   ISZERO(F) returns logical TRUE is F.onefun is identically zero or empty and
%   logical FALSE otherwise.

% Call ISZERO() of the ONEFUN:
if ( isempty(f) )
    out = true;
else
    out = iszero(f.onefun);
end

end
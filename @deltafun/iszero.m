function out = iszero(f)
%ISZERO   True for zero DELTAFUN objects.
%   ISZERO(F) returns logical TRUE if F has a non-zero smooth part of some
%   non trivial delta functions

if ( ~iszero(f.funPart) )
    out = 1;
    return
end

if ( isempty(f.location ) || isempty(f.impulses) )
    out = 1;
    return
end

if ( max(abs(f.impulses(:))) < deltafun.pref.deltafun.deltaTol )
    out = 1;
else
    out = 0;
end

end

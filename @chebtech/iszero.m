function isz = iszero(f)
%ISZERO    True for zero CHEBTECH objects.
%   ISZERO(F) returns logical TRUE is F.values has only zero entries and logical
%   FALSE otherwise.

isz = ~any(f.values, 1);

if ( any(isz) )
    % We need this as any(NaN) = 0, which will pass the test above.
    isz = isz & ~any(isnan(f.values), 1);
    return
end

end
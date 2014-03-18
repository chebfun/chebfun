function out = iszero(f)
%ISZERO    True for zero FOURIERTECH objects.
%   ISZERO(F) returns logical TRUE is F.values has only zero entries and logical
%   FALSE otherwise.

out = ~any(f.values, 1);

if ( any(out) )
    % We need this as any(NaN) = 0, which will pass the test above.
    out = out & ~any(isnan(f.values), 1);
    return
end

end
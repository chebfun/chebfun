function out = iszero(f)
%ISZERO    True for zero CHEBTECH objects.
%   ISZERO(F) returns logical TRUE is F.COEFFS has only zero entries and logical
%   FALSE otherwise.

out = ~any(f.coeffs, 1);

if ( any(out) )
    % We need this as any(NaN) = 0, which will pass the test above.
    out = out & ~any(isnan(f.coeffs), 1);
    return
end

end
function f = abs(f)
%ABS   Absolute value of a CHEBTECH object.
%   ABS(F) returns the absolute value of F, where F is a CHEBTECH object with no
%   roots in F.domain. If ~isempty(roots(F)), then ABS(F) will return garbage
%   with no warning. F may be complex.

if ( isreal(f) || isreal(1i*f) )    
    f.values = abs(f.values);
    f.coeffs = f.chebpoly(f.values);
else
    f = real(sqrt(conj(f).*f));
end

end
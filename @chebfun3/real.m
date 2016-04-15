function f = real(f)
%REAL      Real part of a CHEBFUN3.
%
% See also IMAG.

% Empty check: 
if ( isempty(f) )
    return
end

f = compose(f, @real);

end
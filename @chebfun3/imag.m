function f = imag(f)
%IMAG   Imaginary part of a CHEBFUN3.
%   IMAG(F) returns the imaginary part of a CHEBFUN3.

% Empty check: 
if ( isempty(f) )
    return
end

f = compose(f, @imag); 

end
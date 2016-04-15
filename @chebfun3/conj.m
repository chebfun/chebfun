function f = conj(f)
%CONJ   Complex conjugate of a CHEBFUN3.
%   CONJ(F) returns the complex conjugate of F.  For a complex F, 
%   CONJ(F) = REAL(F) - i*IMAG(F).
%
% See also REAL, IMAG. 

% Check for empty SEPARABLEAPPROX. 
if ( isempty(f) )  
   return
end

% TODO: Write down the formulas for conj, instead of calling the 
% constructor.

f = compose(f, @conj);

end
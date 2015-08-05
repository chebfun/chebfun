function f = conj(f)
%CONJ   Complex conjugate of a SEPARABLEAPPROX.
%   CONJ(F) returns the complex conjugate of F.  For a complex F, CONJ(F) =
%   REAL(F) - i*IMAG(F).
%
% See also REAL, IMAG. 

% Copyright 2015 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% Check for empty SEPARABLEAPPROX. 
if ( isempty( f ) )  
   return
end

% TODO: Write down the formulas for conj, instead of calling the 
% constructor.

f = compose( f, @conj ); 

end

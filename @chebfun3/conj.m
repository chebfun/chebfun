function f = conj(f)
%CONJ   Complex conjugate of a CHEBFUN3.
%   CONJ(F) returns the complex conjugate of F.  For a complex F, 
%   CONJ(F) = REAL(F) - i*IMAG(F).
%
% See also CHEBFUN3/REAL and CHEBFUN3/IMAG.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Check for empty CHEBFUN3.
if ( isempty(f) )  
   return
end

% TODO: Write down the formulas for conj, instead of calling the 
% constructor.

f = compose(f, @conj);

end
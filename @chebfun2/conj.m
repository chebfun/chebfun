function f = conj(f)
%CONJ   Complex conjugate of a CHEBFUN2.
%   CONJ(F) returns the complex conjugate of F.  For a complex F, CONJ(F) =
%   REAL(F) - i*IMAG(F).
%
% See also REAL, IMAG. 

% Copyright 2014 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% Check for empty CHEBFUN2. 
if ( isempty( f ) )  
   return
end

% TODO: Write down the formulas for conj, instead of calling the 
% edconstructor.

op = @(x,y) conj( feval(f, x, y) );  % Resample. 
f = chebfun2( op, f.domain );        % Call constructor. 

end
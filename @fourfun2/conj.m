function f = conj(f)
%CONJ   Complex conjugate of a CHEBFUN2.
%   CONJ(F) returns the complex conjugate of F.  For a complex F, CONJ(F) =
%   REAL(F) - i*IMAG(F).
%
% See also REAL, IMAG. 

% Check for empty CHEBFUN2. 
if ( isempty( f ) )  
   return
end

% TODO: Write down the formulas for conj, instead of calling the 
% edconstructor.

op = @(x,y) conj( feval(f, x, y) );  % Resample. 
f = fourfun2( op, f.domain );        % Call constructor. 

end
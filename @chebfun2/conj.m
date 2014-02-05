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

% TODO: Why is it necessary to resample?

op = @(x,y) conj( feval(f, x, y) );  % Resample. 
f = chebfun2( op, f.domain );        % Call constructor. 

end
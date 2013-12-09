function f = conj(f)
%CONJ Complex conjugate of a chebfun2.
% 
% CONJ(F) returns the complex conjugate of F.  For a complex F, CONJ(F) = 
% REAL(F) - i*IMAG(F). 
%
% See also REAL, IMAG. 

if ( isempty( f ) )  % check for empty chebfun2. 
   return
end

op = @(x,y) conj( feval(f, x, y) );  % Resample. 
f = chebfun2( op, f.domain );               % Call constructor. 

end
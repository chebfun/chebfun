function f = imag( f )
%IMAG   Imaginary part of a CHEBFUN2.
%   IMAG(F) returns the imaginary part of a CHEBFUN2.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Empty check: 
if ( isempty( f ) )
    return
end

op = @(x,y) imag( feval(f, x, y ) );    % Resample.
f = chebfun2( op, f.domain );           % Call constructor. 

end

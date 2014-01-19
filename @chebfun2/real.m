function f = real(f)
%REAL  real part of a chebfun2.
%
% See also IMAG.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.

% Empty check: 
if ( isempty( f ) ) 
    f = chebfun2;
    return
end

op = @(x,y) real( feval( f, x, y ) );  % Resample.
f = chebfun2( op, f.domain );           % Call constructor.

end
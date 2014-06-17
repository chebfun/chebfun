function f = real( f )
%REAL  Real part of a CHEBFUN2.
%
% See also IMAG.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Empty check: 
if ( isempty( f ) ) 
    f = chebfun2;
    return
end

op = @(x,y) real( feval( f, x, y ) );   % Resample.
f = chebfun2( op, f.domain );           % Call constructor.

end

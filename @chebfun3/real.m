function f = real(f)
%REAL      Real part of a CHEBFUN3.
%
% See also IMAG.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Empty check: 
if ( isempty(f) )
    return
end

f = compose(f, @real);

end
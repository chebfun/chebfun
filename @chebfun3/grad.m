function varargout = grad(f)
%GRAD   Gradient of a CHEBFUN3.
%  This command is shorthand for GRADIENT(F).
%
% See also CHEBFUN3/GRADIENT.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Call GRADIENT:
if ( nargout <= 1 )
    out = gradient(f);
    varargout = {out};
else
    [fx, fy, fz] = gradient(f);
    varargout = {fx, fy, fz};
end

end
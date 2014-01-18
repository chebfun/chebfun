function f = abs(f, varargin)
%ABS   Absolute value of a CHEBTECH object.
%   ABS(F) returns the absolute value of F, where F is a CHEBTECH object with no
%   roots in [-1 1]. If ~isempty(roots(F)), then ABS(F) will return garbage
%   with no warning. F may be complex.

%  Copyright 2013 by The University of Oxford and The Chebfun Developers.
%  See http://www.chebfun.org for Chebfun information.

if ( isreal(f) || isreal(1i*f) )    
    f.values = abs(f.values);
    f.coeffs = f.vals2coeffs(f.values);
else
    f = compose(f, @abs, [], varargin{:});
end

end
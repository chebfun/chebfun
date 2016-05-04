function f = abs(f, varargin)
%ABS   Absolute value of a TRIGTECH object.
%   ABS(F) returns the absolute value of F, where F is a TRIGTECH object with no
%   roots in [-1 1]. If ~isempty(roots(F)), then ABS(F) will return garbage
%   with no warning. F may be complex.

%  Copyright 2016 by The University of Oxford and The Chebfun Developers.
%  See http://www.chebfun.org/ for Chebfun information.

if ( isreal(f) || isreal(1i*f) )    
    % Convert values and then compute ABS(). 
    f.values = abs(f.values);
    f.coeffs = f.vals2coeffs(f.values);
    f.isReal = true(1, size(f.coeffs, 2));
else
    f = compose(f, @abs, [], varargin{:});
end

end

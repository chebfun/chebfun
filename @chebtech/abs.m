function f = abs(f, varargin)
%ABS   Absolute value of a CHEBTECH object.
%   ABS(F) returns the absolute value of F, where F is a CHEBTECH object with no
%   roots in [-1 1]. If ~isempty(roots(F)), then ABS(F) will return garbage
%   with no warning. F may be complex.

%  Copyright 2015 by The University of Oxford and The Chebfun Developers.
%  See http://www.chebfun.org/ for Chebfun information.

if ( isreal(f) || isreal(1i*f) )    
    % Convert to values and then compute ABS(). 
    values = f.coeffs2vals(f.coeffs); 
    values = abs(values);
    f.coeffs = f.vals2coeffs(values);
else
    f = compose(f, @abs, [], [], varargin{:});
end

end

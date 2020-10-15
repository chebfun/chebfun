function varargout = coeffs2( f, varargin) 
%COEFFS2   componentwise Chebyshev coefficients of a CHEBFUN2V. 
% 
%   [X, Y] = COEFFS2( F ) returns the bivariate coefficients for each 
%   component of the CHEBFUN2V F in the Chebyshev basis.
%
%   [X, Y] = COEFFS2(F, M, N) returns bivariate coefficients with N Chebyshev 
%   coefficients in the y direction and M Chebyshev coefficients in
%   x direction.

% Copyright 2019 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.


if isempty(f)
    varargout = {}; 
    return
end

F1 = f.components{1};
F2 = f.components{2}; 
if isempty(varargin)
    varargout{1} = coeffs2(F1); 
    varargout{2} = coeffs2(F2); 
else
    varargout{1} = coeffs2(F1, varargin{:}); 
    varargout{2} = coeffs2(F2, varargin{:}); 
end
    

end
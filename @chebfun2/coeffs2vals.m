function varargout = coeffs2vals( U, varargin )
% VAL2COEFFS   Convert matrix of Chebyshev coefficients to values.
%
% V = COEFFS2VALS( C ) given a matrix of bivariate coefficient evaluate the 
% Chebyshev series on a grid size(C).
% 
% [U, S, V] = COEFFS2VALS( U, S, V ) the same as above but keeps everything in low
% rank form. 
% 
% See also VALS2COEFFS.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

tech = chebfunpref().tech(); 
if ( nargin == 1 )
    U = tech.coeffs2vals( flipud(U) ); 
    U = tech.coeffs2vals( (fliplr(U)).' ).'; 
    varargout = {U}; 
    
elseif ( narargin == 3 )
    S = varargin{1}; 
    V = varargin{2}; 
    U = tech.coeffs2vals( flipud(U) ); 
    V = tech.coeffs2vals( (fliplr(V)).' ).'; 
    varargout = {U S V};
    
else
    error('CHEBFUN:CHEBFUN2:coeffs2vals:inputs', ...
        'The number of input arguments should be one or two.');
    
end

end

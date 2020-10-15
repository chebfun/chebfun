function varargout = coeffs2vals( U, varargin )
% VAL2COEFFS    Convert matrix of Chebyshev coefficients to values.
%
% V = COEFFS2VALS( C ) converts a matrix C of bivariate Chebyshev coefficients
% to a matrix of samples V corresponding to values of 
%   sum_i sum_j C(i,j) T_{i-1}(y)T_{j-1}(x) 
% at a tensor Chebyshev grid of size size(C).
% 
% [U, S, V] = COEFFS2VALS( C, D, R ) the same as above but keeps 
% everything in low rank form. Here, C*D*R.' is a 2D matrix of coefficients.
% 
% See also VALS2COEFFS.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( nargin == 1 )
    U = chebtech2.coeffs2vals( U ); 
    U = chebtech2.coeffs2vals( U.' ).'; 
    varargout = {U}; 
    
elseif ( nargin == 3 )
    S = varargin{1}; 
    V = varargin{2}; 
    %U = chebtech2.coeffs2vals( flipud(U) ); 
    %V = chebtech2.coeffs2vals( (fliplr(V)).' ).'; 
    U = chebtech2.coeffs2vals( U ); 
    V = chebtech2.coeffs2vals( V );
    %varargout = {U S V};
    varargout{1} = U; 
    varargout{2} = S; 
    varargout{3} = V; 
    
else
    error('CHEBFUN:CHEBFUN2:coeffs2vals:inputs', ...
        'The number of input arguments should be one or two.');
    
end

end

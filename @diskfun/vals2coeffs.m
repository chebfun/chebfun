function varargout = vals2coeffs( U, varargin )
%VALS2COEFFS    Convert matrix of values to Chebyshev-Fourier coefficients. 
% 
% V = VALS2COEFFS( U ) converts a matrix U of values representing 
% samples of a function from a tensor Chebyshev-Fourier grid 
% to a matrix V of Chebyshev-Fourier coefficients for the corresponding 
% interpolant.
% 
% [U, S, V] = VALS2COEFFS( U, S, V ) the same as above but keeps 
% everything in low rank form.
%
% See also COEFFS2VALS

% Copyright 2019 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.


if ( nargin == 1 )
    U = chebtech2.vals2coeffs( U ); 
    U = trigtech.vals2coeffs( U.' ).'; 
    varargout = {U}; 
elseif ( nargin == 3 )
    S = varargin{1}; 
    V = varargin{2}; 
    U = chebtech2.vals2coeffs( U ); 
    V = trigtech.vals2coeffs( V );     
    varargout = { U S V };
else
    error('CHEBFUN:DISKFUN:vals2coeffs:inputs', ...
        'The number of input arguments should be one or two.');
end

end

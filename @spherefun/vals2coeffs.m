function varargout = vals2coeffs( U, varargin )
%VALS2COEFFS    Convert matrix of values to Fourier-Fourier coefficients. 
% 
% V = VALS2COEFFS( U ) converts a matrix U of values representing 
% samples of a function from a tensor double Fourier sphere grid 
% to a matrix V of 2D Fourier coefficients for the corresponding 
% interpolant.
% 
% [U, S, V] = VALS2COEFFS( U, S, V ) the same as above but keeps 
% everything in low rank form. Here, U*S*V.' is a sample of 
% values on the double Fourier sphere grid.
% 
% See also COEFFS2VALS

% Copyright 2019 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.


if ( nargin == 1 )
    U = trigtech.vals2coeffs( U ); 
    U = trigtech.vals2coeffs( U.' ).'; 
    varargout = {U}; 
elseif ( nargin == 3 )
    S = varargin{1}; 
    V = varargin{2}; 
    U = trigtech.vals2coeffs( U ); 
    V = trigtech.vals2coeffs( V );     
    varargout = {U S V};
else
    error('CHEBFUN:SPHEREFUN:vals2coeffs:inputs', ...
        'The number of input arguments should be one or two.');
end

end

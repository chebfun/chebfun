function varargout = coeffs2vals( U, varargin )
% COEFFS2VALS    Convert matrix of Chebyshev-Fourier coefficients to values.
%
% V = COEFFS2VALS( U ) converts a matrix U of Chebyshev-Fourier
% coefficients to a matrix of samples V corresponding to values sampled
% on a doubled up tensor polar grid of Chebyshev points in 
% the radial direction from [-1, 1], and equispaced points in the
% angular direction from [-pi, pi). 
%
% [U, S, V] = COEFFS2VALS( C, D, R ) the same as above but keeps everything in
% low rank form. Here, C*D*R.' is a 2D coeff matrix.
%

% See also VALS2COEFFS.

% Copyright 2019 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

    
if ( nargin == 1 )
    
    % Get size: 
    [m, n] = size( U ); 

    
    U = chebtech2.coeffs2vals( U ); 
    U = trigtech.coeffs2vals( U.' ).'; 
    U = real(U); %assume function is real-valued
    varargout = {U}; 
    
elseif ( nargin == 3 )
    S = varargin{1}; 
    V = varargin{2}; 
    U = chebtech2.coeffs2vals( U ); 
    V = trigtech.coeffs2vals(V); 
    varargout = {U S V};
    
else
    error('CHEBFUN:DISKFUN:coeffs2vals:inputs', ...
        'The number of input arguments should be one or three.');
    
end

   
end

function varargout = coeffs2vals( U, varargin )

if ( nargin == 1 )
    U = chebtech2.coeffs2vals( U ); 
    U = chebtech2.coeffs2vals( U' ); 
    varargout = {U}; 
elseif ( narargin == 3 )
    S = varargin{1}; 
    V = varargin{2}; 
    U = chebtech2.coeffs2vals( U ); 
    V = chebtech2.coeffs2vals( V ); 
    varargout = {U S V};
else
    error('CHEBFUN2:COEFFS2VALS:inputs', 'The number of input arguments should be one or two');
end
function varargout = vals2coeffs( U, varargin )

if ( nargin == 1 )
    U = chebtech2.vals2coeffs( U ); 
    U = chebtech2.vals2coeffs( U' ); 
    varargout = {U}; 
elseif ( narargin == 3 )
    S = varargin{1}; 
    V = varargin{2}; 
    U = chebtech2.vals2coeffs( U ); 
    V = chebtech2.vals2coeffs( V ); 
    varargout = {U S V};
else
    error('CHEBFUN2:VALS2COEFFS:inputs', 'The number of input arguments should be one or two');
end
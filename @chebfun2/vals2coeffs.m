function varargout = vals2coeffs( U, varargin )
% VALS2COEFFS   Convert matrix of values to Chebyshev coefficients. 
% 
% See also COEFFS2VALS

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.

if ( nargin == 1 )
    U = chebtech2.vals2coeffs( U ); 
    U = chebtech2.vals2coeffs( U.' ).'; 
    varargout = {U}; 
elseif ( narargin == 3 )
    S = varargin{1}; 
    V = varargin{2}; 
    U = chebtech2.vals2coeffs( U ); 
    V = chebtech2.vals2coeffs( V.' ).'; 
    varargout = {U S V};
else
    error('CHEBFUN2:VALS2COEFFS:inputs', 'The number of input arguments should be one or two');
end
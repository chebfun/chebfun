function varargout = coeffs2vals( U, varargin )
% VAL2COEFFS   Convert matrix of Chebyshev coefficients to values. 
% 
% See also VALS2COEFFS.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.

% TODO: Document.

if ( nargin == 1 )
    U = chebtech2.coeffs2vals( U ); 
    U = chebtech2.coeffs2vals( U' ).'; % TODO: U.' ?
    varargout = {U}; 
    
elseif ( narargin == 3 )
    S = varargin{1}; 
    V = varargin{2}; 
    U = chebtech2.coeffs2vals( U ); 
    V = chebtech2.coeffs2vals( V.' ).'; 
    varargout = {U S V};
    
else
    error('CHEBFUN2:COEFFS2VALS:inputs', ...
        'The number of input arguments should be one or two.');
    
end

end
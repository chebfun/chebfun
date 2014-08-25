function varargout = vals2coeffs( U, S, V )
%VALS2COEFFS   Convert matrix of values to Chebyshev coefficients. 
% 
% V = VALS2COEFFS( C ) given a matrix C of va.ues on a tensor grid, this returns
% the corresponding bivaraite Chebyshev coefficients, V, which is a matrix of 
% size(C).
% 
% [U, S, V] = VALS2COEFFS( U, S, V ) the same as above but keeps everything in low
% rank form. 
% 
% See also COEFFS2VALS

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

tech = chebfunpref().tech();

if ( nargin == 1 )
    U = tech.vals2coeffs( U ); 
    U = tech.vals2coeffs( flipud(U).' ).'; 
    U = fliplr(U);
    varargout = {U}; 
elseif ( narargin == 3 )
    U = tech.vals2coeffs( U ); 
    V = tech.vals2coeffs( V.' ).'; 
    U = fliplr(U);
    V = fliplr(V);
    varargout = {U S V};
else
    error('CHEBFUN:CHEBFUN2:vals2coeffs:inputs', ...
        'The number of input arguments should be one or two.');
end

end

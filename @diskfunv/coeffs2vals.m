function [X, Y] = coeffs2vals(U,V)
% COEFFS2VALS  componentwise conversion of matrix of Chebyshev-Fourier 
% coefficients to values.
%
%  X,Y = COEFFS2VALS( U, V ) converts matrices U and V of Chebyshev-Fourier
% coefficients to matrices of samples X and Y corresponding to values sampled
% on a doubled up tensor polar grid of Chebyshev points in 
% the radial direction from [-1, 1], and equispaced points in the
% angular direction from [-pi, pi). 
%
% See also COEFFS2, VALS2COEFFS

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information. 

%%

X = diskfun.coeffs2vals(U); 
Y = diskfun.coeffs2vals(V); 
end




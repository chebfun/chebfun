function F = coeffs2diskfun(X)
%COEFFS2DISKFUN   Convert a matrix of Cheb-Fourier coefficients to a diskfun. 
% 
%  F = coeffs2diskfun(X) returns a diskfun object F that has a
%  Fourier--Chebyshev matrix of coefficients X.  This is useful for
%  computing quantities on the disk with the function F.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Get size: 
[m, n] = size( X ); 

% If n is odd, then make it even. 
if ( mod(n, 2) == 1 ) 
    X = [ zeros(m, 1) X ]; 
end

% Convert to values on the grid: 
vals = trigtech.coeffs2vals(chebtech2.coeffs2vals(X).').'; 

% Assume that the function is real-valued:
vals = real( vals );

% Restrict to the region of interest: 
vals = vals(floor(m/2)+1:m, :);
 

% Finally, make a diskfun object out of the values: 
F = diskfun( vals );   % constructor

end
function F = coeffs2diskfun(X)
%COEFFS2DISKFUN   Convert a matrix of Chebyshev-Fourier coefficients to a 
%                 diskfun. 
% 
%   F = coeffs2diskfun( X ) returns a diskfun object F that has a
%   Chebyshev-Fourier matrix of coefficients X.  This is useful for
%   computing quantities on the disk with the function F.
% 
% See also DISKFUN/COEFFS2

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% To get the correct doubled-up grid, we can only accept
% certain coeff matrix sizes: 

% Get size: 
[m, n] = size( X ); 

% If n is odd, then make it even. 
if ( mod(n, 2) == 1 ) 
    X = [ zeros(m, 1) X ]; 
end

%If m is even, we need it to be odd: 
if (mod(m, 2) == 0)
    m=m+1;
    X = chebtech2.alias(X, m); 
end

% Convert to values on the grid.
vals = trigtech.coeffs2vals(chebtech2.coeffs2vals(X).').'; 

% Assume that the function is real-valued.
vals = real( vals );

% Restrict to the region of interest.
vals = vals(floor(m/2)+1:m, :);
 
% Finally, make a diskfun object out of the values.
F = diskfun( vals );

end
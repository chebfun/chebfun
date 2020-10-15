function F = coeffs2spherefun( X )
%COEFFS2SPHEREFUN      Make a SPHEREFUN object from a matrix of Fourier
%   coefficients. 
%
%  F = coeffs2spherefun(X) returns a spherefun object F that has a
%  Fourier--Fourier matrix of coefficients X.  This is useful for
%  computing quantities on the sphere with the function F.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Get size: 
[m, n] = size( X ); 

% If n is odd, then make it even. 
if ( mod(n, 2) == 1 ) 
    X = [ zeros(m, 1) X ];
    n = n+1;
end

if ( mod(m, 2) == 1 ) 
    X = [ zeros(1, n); X  ];
    m = m + 1; 
end

% Convert to values on the grid: 
VALS = trigtech.coeffs2vals(trigtech.coeffs2vals(X).').'; 
% Restrict to the region of interest: 
VALS = VALS([floor(m/2)+1:m 1], :);

% Only real values are supported.
VALS = real(VALS); 

% Finally, make a spherefun object out of the values: 
F = spherefun(VALS); 

end
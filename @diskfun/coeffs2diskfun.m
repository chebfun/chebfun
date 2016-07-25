function F = coeffs2diskfun(X)
%COEFFS2DISKFUN      Make a DISKFUN object from a matrix of
%Fourier-Chebyshev coefficients. 
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

%if ( mod(m, 2) == 1 ) %m can be/should be odd.
 %   X = [ zeros(1, n+1); X  ];
  %  m = m + 1; 
%end

% Convert to values on the grid: 
VALS = trigtech.coeffs2vals(chebtech2.coeffs2vals(X).').'; 
VALS = real(VALS);
% Restrict to the region of interest: 
VALS = VALS(floor(m/2)+1:m, :);



% Finally, make a diskfun object out of the values: 
F = diskfun(VALS); 

end
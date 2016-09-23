function v = std2(f)
%STD2   Standard deviation of a SEPARABLEAPPROX.
%   V = STD2(F) computes the standard deviation of a SEPARABLEAPPROX, i.e., 
%
%     STD2(F)^2 = 1/A*sum2(|f(x,y) - m|^2)
%
%   where A is the area of the domain of F.
%
% See also MEAN, MEAN2, STD.

% Copyright 2016 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

h = f - mean2( f ); 
v = sqrt( mean2 (h .* conj( h ) ) );

end

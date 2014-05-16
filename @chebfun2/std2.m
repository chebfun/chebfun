function v = std2(f)
%STD2   Standard deviation of a CHEBFUN2.
%   V = STD2(F) computes the standard deviation of a CHEBFUN2, i.e., 
%                                 d  b
%                                 /  /
%     STD2(F)^2 = 1/(b-a)/(d-c)  |  |  |f(x,y) - m|^2 dx dy
%                                /  /
%                               c  a 
%   where the domain of F is [a,b] x [c,d], and m = mean2(F). 
%
% See also MEAN, MEAN2, STD.

% Copyright 2014 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

h = f - mean2( f ); 
v = sqrt( mean2 (h .* conj( h ) ) );

end
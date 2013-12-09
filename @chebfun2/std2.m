function v = std2(f)
% STD2 Standard deviation of a chebfun2. 
%
% v = STD2(F) computes the standard deviation of a chebfun2, i.e., 
% 
%                             d  b
%                             /  /
%  STD2(F)^2 = 1/(b-a)/(d-c)  |  |  |f(x,y) - m|^2 dx dy
%                            /  /
%                           c  a 
% 
% where the domain of F is [a,b] x [c,d], and m = mean2(F). 
%
% See also MEAN, MEAN2, STD.

h = f - mean2( f ); 
v = sqrt( mean2 (h .* conj( h ) ) );

end
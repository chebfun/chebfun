function v = mean2(f)
%MEAN2 mean of a chebfun2
%
% V = MEAN2(F) returns the mean of a chebfun2: 
% 
%                       d  b
%                      /  /   
%  V = 1/(d-c)/(b-a)   |  |   f(x,y) dx dy 
%                      /  /
%                     c  a
% 
% where the domain of f is [a,b] x [c,d]. 
%
% See also MEAN, STD2.

dom = f.domain; 
width = diff( dom(1:2) ); 
height = diff( dom(3:4) ); 
area = width * height; 
v = sum2( f ) / area; 

end
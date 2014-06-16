function v = mean2( f )
%MEAN2   Mean of a CHEBFUN2
%   V = MEAN2(F) returns the mean of a CHEBFUN: 
% 
%                        d  b
%                       /  /   
%   V = 1/(d-c)/(b-a)   |  |   f(x,y) dx dy 
%                       /  /
%                      c  a
% 
% 	where the domain of F is [a,b] x [c,d]. 
%
% See also MEAN, STD2.

% Copyright 2014 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% Empty check:
if ( isempty( f ) ) 
    return
end 

% Apply the formula: 
dom = f.domain; 
width = diff( dom(1:2) ); 
height = diff( dom(3:4) );   
area = width * height; 
v = sum2( f ) / area;  

end
function v = mean2( f )
%MEAN2   Mean of a SEPARABLEAPPROX.
%   V = MEAN2(F) returns the mean of a SEPARABLEAPPROX: 
% 
%   V = 1/A*integral2( f )
% 
% 	where the A is the area of the domain of F. 
%
% See also MEAN, STD2.

% Copyright 2015 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% Empty check:
if ( isempty( f ) ) 
    return
end 

% Apply the formula: 
v = sum2( f ) / domainarea( f );  

end

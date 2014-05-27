function map = linear2D(corners)
% MAP = LINEAR2D(CORNERS) creates a map structure for funs2
% it returns a structure that defines a linear map.  
%
% MAP.FOR is a function that maps [-1,1]^2 to 
% [corners(1) corners(2)] x [corners(3) corners(4)].
% MAP.INV is the inverse map. 
% MAP.ID is a string that identifies the map. 

% Extract out input. 
a = corners(1); b = corners(2); c=corners(3); d=corners(4);

% linear2D map:
map = struct('for', @(s,t) [b*(s+1)/2+a*(1-s)/2,d*(t+1)/2+c*(1-t)/2],...
             'inv', @(x,y) [(x-a)/(b-a)-(b-x)/(b-a), (y-c)/(d-c)-(d-y)/(d-c)],...
             'derx',@(y) (b-a)/2 + 0*y,'dery',@(y) (d-c)/2 + 0*y,'name','linear','par',[a b c d]) ;

end
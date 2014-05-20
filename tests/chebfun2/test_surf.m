function pass = test_surf( pref ) 
% Test surf

if ( nargin == 0) 
    pref = chebfunpref; 
end


f = chebfun2(@(x,y) cos(x.*y)); 
x = chebfun2(@(x,y) x); 
y = chebfun2(@(x,y) y); 

pass = 1; 
try 
   surf( f )
   surf( f , f)
   surf( f, f, f )
   surf(f, 'numpts', 10)
   surf(x, y, f)
   surf(x, y, f, 'edgecolor', 'r')
catch ME 
    pass = 0; 
end
close all

end
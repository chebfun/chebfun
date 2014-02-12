function pass = test_chebfun2_contour( pref ) 
% Test contour

if ( nargin == 0) 
    pref = chebpref; 
end

tol = 1000*pref.cheb2Prefs.eps; 
j = 1; 

f = chebfun2(@(x,y) cos(x.*y)); 
x = -1:.1:1;
[xx, yy] = meshgrid(x); 

pass = 1; 
try 
   contour( f )
   contour( f , 5)
   contour( f, [0 0])
   contour(f, 'numpts', 100)
   contour(f, 'pivots', 'r.-')
   contour(xx, yy, f)
   contourf( f )
   contourf( f , 5)
   contourf( f, [0 0])
   contourf(f, 'numpts', 100)
   contourf(f, 'pivots', '.')
   contourf(xx, yy, f)
catch ME 
    pass = 0; 
end

close all
end
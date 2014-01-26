function pass = test_chebfun2_surf( pref ) 
% Test surf

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
   surf( f )
   surf( f , f)
   surf( f, f, f )
   surf(f, 'numpts', 100)
   surf(xx, yy, f)
   surf(xx, yy, f, 'edgecolor', 'r')
   surfc( f )
   surfc( f , f)
   surfc( f, f, f )
   surfc(f, 'numpts', 100)
   surfc(xx, yy, f)
   surfc(xx, yy, f, 'edgecolor', 'r')
catch ME 
    pass = 0; 
end
end
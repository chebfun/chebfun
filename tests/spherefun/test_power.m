function pass = test_power( ) 
% Test power in SPHEREFUN 

tol = 1e3*chebfunpref().techPrefs.eps;

f = spherefun(@(x,y,z) z );
g = spherefun(@(x,y,z) z.^2 );
pass(1) = norm( f.^2 - g ) < tol; 

end 
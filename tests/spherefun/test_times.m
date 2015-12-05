function pass = test_times( ) 
% Test times in SPHEREFUN 

tol = 1e3*chebfunpref().techPrefs.eps;

f = spherefun(@(x,y,z) sin(x.*y.*z)); 
pass(1) = norm( f.*f - f.^2 ) < tol; 

end 
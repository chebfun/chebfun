function pass = test_mean( ) 
% Test MEAN()

tol = 1000*chebfunpref().cheb2Prefs.chebfun2eps;

g = diskfun();
pass(1) = isempty(mean(g)); 

%%
g = diskfun(@(x,y) 0*x+1); 
pass(2) = ( norm(mean(g)-1) < tol );
pass(3) = ( norm(mean(g,2)-1)<tol );

g = diskfun(@(x,y) exp(-10*x.^2-10*(y-.3).^2)-10);
pass(4) = ( norm(mean(g)-sum(g)) < tol ); 
pass(5) = ( norm(mean(g,2) - sum(g, 2)/(2*pi)) < tol); 
end
function pass = test_sum( ) 
% Test SUM()

tol = 1000*chebfunpref().cheb2Prefs.chebfun2eps;

g = diskfun();
pass(1) = isempty(sum(g)); 
pass(2) = isempty(sum(g, 2));

%%
g = diskfun(@(x,y) 0*x+1); 
pass(3) = ( norm(sum(g)-1/2) < tol ); %test sum over angular vars 
pass(4) = ( norm(sum(g,2)-2*pi) < tol ); %test sum over radial vars

end
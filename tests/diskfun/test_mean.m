function pass = test_mean( ) 
% Test MEAN()

tol = 1000*chebfunpref().cheb2Prefs.chebfun2eps;

g = diskfun();
pass(1) = isempty(mean(g)); 

%%
g = diskfun(@(x,y) 0*x+1); 
pass(2) = ( norm(mean(g)-.5) < tol ); %this direction includes measure on disk
pass(3) = ( norm(mean(g,2)-1)<tol ); 

g = diskfun(@(t,r) r.^2.*sin(2*t), 'polar');
pass(4) = ( norm(mean(g)-chebfun(@(t) .25*sin(2*t), [-pi pi], 'trig') ) < tol );
pass(5) = ( norm(mean(g,2)) < tol); 

g = diskfun( @(x,y) exp(-10*x.^2-10*y.^2));
pass(6) = ( abs(sum(sum(g)) - sum2(g)) < tol ); 
end
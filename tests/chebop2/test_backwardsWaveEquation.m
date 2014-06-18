function pass = test_backwardsWaveEquation( prefs )
% Check that the backwards wave equation is working. 
% Alex Townsend, August 2013. 

if ( nargin < 1 ) 
    prefs = chebfunpref(); 
end 
tol = 100*prefs.techPrefs.eps; 

%%

d = [-pi pi 0 1]; 
exact = chebfun2(@(x,t) sin(x+t), d); 
N = chebop2(@(u) diff(u,2,1) - diff(u,2,2), d);
N.lbc = @(t) sin(-pi+t);
N.rbc = @(t) sin(pi+t);
N.ubc = @(x,u) [u-sin(x+1) ; diff(u)-cos(x+1)];
u = N \ 0;
pass(1) = ( norm(u - exact) < 10*tol ) ;

end
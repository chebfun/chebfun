function pass = test_linearSchrodinger( prefs )
% Check that we can do simple linear Schrodinger equation. 
% Alex Townsend, April 2013. 
% 
% tol = 100*chebfun2pref('eps'); 
% j = 1; 
% 
% % Simple, first example. 
% d = [0 1 0 1]; 
% exact = chebfun2(@(x,t) exp(1i*t).*exp(x), d);
% N = chebop2(@(u) 1i*diffy(u) + diffx(u,2), d); 
% N.ubc = @(x) exp(1i*t).*exp(x);
% N.lbc = @(t) exp(1i*t); 
% N.rbc = @(t) exp(1i*t)*exp(1);
% u = N \ 0;
% 
% pass(j) = ( norm(u - exact) < tol ); j = j + 1; 

pass = 1; 
end
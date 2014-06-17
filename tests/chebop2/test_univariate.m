function pass = test_univariate( prefs )
% Rank-1 PDEs can be solved by univariate spectral methods.  Do the results
% match. 
% Alex Townsend, March 2013. 

if ( nargin < 1 ) 
    prefs = chebfunpref(); 
end 
tol = 1000*prefs.techPrefs.eps; 

% Simple example in y-variable.
N = chebop2(@(u) diff(u,2,1) + u); N.dbc = 1; N.ubc = 1; u = N \ 0; 
L = chebop(@(u) diff(u,2) + u); L.lbc = 1; L.rbc = 1; v = L \ 0; 

pass(1) = (length(u) == 1);   
pass(2) = (norm(u(0,:) - v) < tol); 


% Simple example in x-variable.
N = chebop2(@(u) diff(u,2,2) + u); N.lbc = 1; N.rbc = 1; u = N \ 0; 
L = chebop(@(u) diff(u,2) + u); L.lbc = 1; L.rbc = 1; v = L \ 0; 

pass(3) = (length(u) == 1);  
pass(4) = (norm(u(:,0) - v.') < tol);  


% 
% % Simple example in x-variable.
% N = chebop2(@(u) diff(diff(u,1,1),1,2) + diff(u,1,1)); 
% N.lbc = 1; N.dbc = 1; u = N \ 0; 
% 
% L = chebop(@(u) diff(u,2) + u); L.lbc = 1; L.rbc = 1; v = L \ 0; 
% 
% pass(j) = (length(u) == 1);  j = j+1; 
% pass(j) = (norm(u(:,0) - v) < tol); j = j+1; 
function pass = test_techs( pref )
% Check that Chebfun2 can work with the different techs. 

if ( nargin < 1 )
    pref = chebfunpref(); 
end

tol = 100*pref.eps;

%%
pref.tech = @chebtech1;
F = @(x,y) cos(x.*y);
f = chebfun2(F, [-1, 1, -1, 1], pref); 
pass(1) = norm(feval(f,0,0) - F(0,0)) < tol;

%%
pref.tech = @chebtech2;
F = @(x,y) cos(x.*y);
f = chebfun2(F, [-1, 1, -1, 1], pref); 
pass(2) = norm(feval(f,0,0) - F(0,0)) < tol;

%%
pref.tech = @trigtech;
F = @(x,y) cos(pi*y).*sin(pi*x);
f = chebfun2(F, [-1, 1, -1, 1], pref); 
pass(3) = norm(feval(f,0,0) - F(0,0)) < 100*tol;

end

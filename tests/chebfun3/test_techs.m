function pass = test_techs(pref)
% Check that Chebfun3 can work with the different techs. 

if ( nargin < 1 )
    pref = chebfunpref(); 
end
tol = 100*pref.cheb3Prefs.chebfun3eps;

%%
pref.tech = @chebtech1;
ff = @(x,y,z) cos(x.*y.*z);
dom = [-1, 1, -1, 1 -1 1];
f = chebfun3(ff, dom, pref);
pass(1) = isa(f.cols.funs{1}.onefun, 'chebtech1');
pass(2) = abs(feval(f, 0, 0, 0) - ff(0, 0, 0)) < tol;

%%
pref.tech = @chebtech2;
f = chebfun3(ff, dom, pref); 
pass(3) = isa(f.cols.funs{1}.onefun, 'chebtech2');
pass(4) = abs(feval(f, 0, 0, 0) - ff(0, 0, 0)) < tol;

%%
pref.tech = @trigtech;
ff = @(x,y,z) cos(pi*z).*sin(pi*y).*sin(pi*x);
f = chebfun3(ff, dom, pref);
pass(5) = isa(f.cols.funs{1}.onefun, 'trigtech');
pass(6) = abs(feval(f, 0, 0, 0) - ff(0, 0, 0)) < 100*tol;

end
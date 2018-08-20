function pass = test_spherefun2ballfun( ) 
% Test with function cos(cos(lam)*sin(th))
S = [21,23,25];
f = spherefun(@(lam,th)cos(cos(lam).*sin(th)));
exact = ballfun(@(r,lam,th)cos(cos(lam).*sin(th)),S);
g = ballfun.spherefun2ballfun(f,S);
pass(1) = isequal(g,exact);
end

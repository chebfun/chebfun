function pass = test_spherefun2ballfun( ) 
% Test with function cos(cos(lam)*sin(th))
f = spherefun(@(lam,th)cos(cos(lam).*sin(th)));
[p,n] = size(coeffs2(f));
exact = ballfun(@(r,lam,th)cos(cos(lam).*sin(th)),[2,n,p]);
g = ballfun.spherefun2ballfun(f);
pass(1) = isequal(g,exact);
end

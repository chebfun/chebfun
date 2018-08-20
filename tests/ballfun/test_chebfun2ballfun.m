function pass = test_chebfun2ballfun( ) 
% Test with function sin(r)*sqrt(r)
f = @(r)sin(r)*sqrt(r);
S = [10,15,42];
g = ballfun.chebfun2ballfun(chebfun(f,S(1)),S);
h = ballfun(@(r,lam,th)sin(r).*sqrt(r),S);

pass(1) = isequal(g,h);
end

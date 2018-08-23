function pass = test_chebfun2ballfun( ) 

% Test with function sin(r)
f = chebfun(@(r)sin(r));
n = length(f);
g = ballfun.chebfun2ballfun(f);
h = ballfun(@(r,lam,th)sin(r));
pass(1) = isequal(g,h);
end

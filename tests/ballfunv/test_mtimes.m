function pass = test_mtimes( ) 
% Test with function coeffs 1
f = ballfun(ones(20,21,22));
F = ballfunv(f,f,f);
g = ballfun(2*ones(20,21,22));
G = ballfunv(g,g,g);
pass(1) = (isequal(2*F,G) && isequal(F*2,G));
end

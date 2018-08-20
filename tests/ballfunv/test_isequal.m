function pass = test_isequal( ) 
% Test with function 1
f = ballfun(ones(20,21,22));
F = ballfunv(f,f,f);
G = F+F-F;

pass(1) = (isequal(F,G) && isequal(G,F));
end

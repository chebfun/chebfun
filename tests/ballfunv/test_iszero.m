function pass = test_iszero( ) 
% Test with vector 1-1
f = ballfun(ones(20,21,22));
F = ballfunv(f,f,f);
G = F-F;
pass(1) = iszero(G);
end

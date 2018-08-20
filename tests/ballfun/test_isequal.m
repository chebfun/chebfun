function pass = test_isequal( ) 
% Test with function 1
f = ballfun(ones(20,21,22));
g = f+f-f;

pass(1) = (isequal(f,g) && isequal(g,f));
end

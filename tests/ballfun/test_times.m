function pass = test_times( ) 
% Test with function 1, 2*3=6
f = ballfun(2*ones(20,21,22),"vals");
g = ballfun(3*ones(20,21,22),"vals");
h = ballfun(6*ones(20,21,22),"vals");

pass(1) = (isequal(f*g,h) && isequal(g*f,h));
end

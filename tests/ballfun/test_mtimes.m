function pass = test_mtimes( ) 
% Test with function coeffs 1
f = ballfun(ones(20,21,22));
g = ballfun(2*ones(20,21,22));

pass(1) = (isequal(2*f,g) && isequal(f*2,g));
end

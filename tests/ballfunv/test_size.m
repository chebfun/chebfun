function pass = test_size( ) 
% Test with function coeffs 1

f = ballfun(ones(20,21,22));
F = ballfunv(f,f,f);
pass(1) = (min(size(F)==[20,21,22]) == 1);

f = ballfun(ones(1000,1,1));
F = ballfunv(f,f,f);
pass(2) = (min(size(F)==[1000,1,1]) == 1);

if (nargout > 0)
    pass = all(pass(:));
end
end

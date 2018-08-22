function pass = test_solHarm()
% Test the equality for P^m_l, 0 <= l <= n and -m <= l <= m
eps = 1e-10;

n = 10;
Max_difference = 0;

for l = 0:n
    for m = 0:l
        Y = extract_spherefun(ballfun.solHarm(l,m))*sqrt(2*l+1)*sqrt(1+(m>0));
        Z = spherefun.sphharm(l,m);
        Max_difference = max(norm(Y-Z),Max_difference);
    end
end
pass = Max_difference < eps;
end

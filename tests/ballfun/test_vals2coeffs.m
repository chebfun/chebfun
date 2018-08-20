function pass = test_vals2coeffs( ) 
% Example 1
X = ones(21,22,23);
V = ballfun.vals2coeffs(X);
V2 = zeros(21,22,23);
V2(1,12,12) = 1;
pass(1) = (max(max(max(abs(V-V2))))<1e-15);

% Example 2
S = [100,150,200];
m = S(1); n = S(2); p = S(3);
X = ones(S);
V = ballfun.vals2coeffs(X);
V2 = zeros(S);
V2(1,floor(n/2)+1,floor(p/2)+1) = 1;
pass(2) = (max(max(max(abs(V-V2))))<1e-15);

if (nargout > 0)
    pass = all(pass(:));
end
end

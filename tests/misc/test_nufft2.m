function pass = test_nufft2( pref ) 

if ( nargin == 0 ) 
    pref = chebfunpref();
end 

rng(0);
M = 100; N = 103; 
m = 105; n = 101; 
C = rand(m, n) + 1i*rand(m,n);

x = rand(M, N); 
y = rand(M, N);

% Fast transform: 
F = chebfun.nufft2( C, x, y ); 

% direct 
Y = zeros(M, N); 
for j = 1:M
    for k = 1:N 
        Y(j,k) = exp(-2*pi*1i*y(j,k)*(0:m-1))*C*exp(-2*pi*1i*x(j,k)*(0:n-1)');
    end
end 

pass(1) = norm( F - Y )/norm(C) < 10*M*N*eps;

end
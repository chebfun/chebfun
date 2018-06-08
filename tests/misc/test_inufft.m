function pass = test_inufft( pref )

if ( nargin == 0 ) 
    pref = chebfunpref();
end

% Choose a tolerance:
tol = eps;
rng(0)

% Test INUFFT-I 
count = 1;
for N = [10.^(0:3) 4000]
    omega = (0:N-1)'; 
    omega = omega + .4*rand(N,1);
    F = exp(-2*pi*1i*(0:N-1)'/N*omega.');
    c = rand(N,1) + 1i*rand(N,1);
    exact = F \ c; 
    fast = chebfun.inufft( c, omega, 1 );
    pass(count) = norm( exact - fast, inf ) < 300*N*tol*norm(c,1);
    count = count + 1;
end

% Check the plan works.
[fast,plan] = chebfun.inufft( c, omega, 1 );
pass(count) = norm( fast - plan(c) ) == 0;
count = count + 1;

% Test INUFFT-II:
for N = [10.^(0:3) 4000]
    x = linspace(0,1,N+1)'; x(end) = [];
    x = x + .2*rand(N,1)/N;
    F = exp(-2*pi*1i*x*(0:N-1)); 
    c = rand(N,1) + 1i*rand(N,1);
    exact = F \ c; 
    fast = chebfun.inufft( c, x, 2 );
    pass(count) = norm( exact - fast, inf ) < 300*N*tol*norm(c,1);
    count = count + 1; 
end

% Check the plan works.
[fast,plan] = chebfun.inufft( c, x, 2 );
pass(count) = norm( fast - plan(c) ) == 0;
count = count + 1;

% Check two input parameters gives INUFFT-II
N = 100;
x = linspace(0,1,N+1)'; x(end) = [];
x = x + .2*rand(N,1)/N;
F = exp(-2*pi*1i*x*(0:N-1)); 
c = rand(N,1) + 1i*rand(N,1);
exact = F \ c; 
fast = chebfun.inufft( c, x );
pass(count) = norm( exact - fast, inf ) < 300*N*tol*norm(c,1);
count = count + 1;

end
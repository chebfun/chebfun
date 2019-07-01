function pass = test_nufft( pref )

if ( nargin == 0 ) 
    pref = chebfunpref();
end

% Choose a tolerance:
tol = eps;
rng(0)

% Test NUFFT-I 
count = 1; 
for N = 10.^(0:4)
    omega = (0:N-1)'; 
    omega = omega + 1.2*rand(N,1)/N;
    c = rand(N,1) + 1i*rand(N,1);
    exact = nudft1( c, omega ); 
    fast = chebfun.nufft( c, omega, 1 );
    pass(count) = norm( exact - fast, inf ) < 300*N*tol*norm(c,1);
    count = count + 1;
end

% Test on random inputs: 
N = 50; 
omega = N*randn(1,N);
c = rand(N,1) + 1i*rand(N,1);
F = exp(-2*pi*1i*((0:N-1)/N).'*omega);
f = chebfun.nufft(c, omega.', 1);
pass(count) = ( norm( f - F*c ) < 300*N*tol*norm(c,1) );
count = count + 1; 

% Test NUFFT-II:
for N = 10.^(0:4)
    x = linspace(0,1,N+1)'; x(end) = [];
    x = x + 1.2/N;
    c = rand(N,1) + 1i*rand(N,1);
    exact = nudft2( c, x ); 
    fast = chebfun.nufft( c, x );
    pass(count) = norm( exact - fast, inf ) < 300*N*tol*norm(c,1);
    count = count + 1; 
end

% Test on random inputs: 
N = 50; 
x = 100*rand(N,1);
c = rand(N,1) + 1i*rand(N,1);
F = exp(-2*pi*1i*x*(0:N-1));
f = chebfun.nufft(c, x);
pass(count) = ( norm( f - F*c ) < 300*N*tol*norm(c,1) );
count = count + 1; 

% Test NUFFT-III
for N = 10.^(0:4)
    omega = (0:N-1)'; 
    omega = omega + 1.2*rand(N,1)/N;
    x = linspace(0,1,N+1)'; x(end) = [];
    x = x + 1.2/N;
    c = rand(N,1) + 1i*rand(N,1);
    exact = nudft3( c, x, omega ); 
    fast = chebfun.nufft( c, x, omega, 3 );
    pass(count) = norm( exact - fast, inf ) < 10000*N*tol*norm(c,1);
    count = count + 1;
end

% Test NUFFT-II on nonsquare inputs: 
n = 101; 
m = 3; 
x = linspace(0,1,m+1)'; x(end) = [];
x = x + 1.01/m;
c = rand(n,1); 
exact = nudft2( c, x );
fast = chebfun.nufft( c, x );
pass(count) = norm( exact - fast, inf ) < 100*tol*norm(c,1);
count = count + 1;

end

function f = nudft1( c, omega) 

f = zeros(size(omega,1),1);
for j = 1:numel(f)
    f(j) = exp(-2*pi*1i*(j-1)/size(omega,1)*omega.')*c;
end
end

function f = nudft2( c, x) 

f = zeros(size(x,1),1);
omega = 0:size(c,1)-1;
for j = 1:numel(f)
    f(j) = exp(-2*pi*1i*x(j)*omega)*c;
end
end

function f = nudft3( c, x, omega) 

f = zeros(size(omega,1),1);
for j = 1:numel(f)
    f(j) = exp(-2*pi*1i*x(j)*omega.')*c;
end
end
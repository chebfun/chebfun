function pass = test_sample( ) 
% Test diskfun sample() command 

tol = 100*chebfunpref().cheb2Prefs.chebfun2eps;

% Function to test
f = diskfun(@(x,y) sin(pi*x.*y));

% Ensure the matrix of sampled values is correct.
[m,n] = length(f);
[nn,mm] = size(sample(f));
pass(1) = (m == mm) && (n == nn);

% Sample on fixed grids of various sizes to make sure the right size output
% is given.
m = 120; 
n = 121;
[nn,mm] = size(sample(f, m, n));
pass(2) = (m == mm) && (n == nn);

m = 121; 
n = 120;
[nn,mm] = size(sample(f, m, n));
pass(3) = (m == mm) && (n == nn);

% Check samples are correct.
% m even and n odd
m = 30; 
n = 2*20+1;
cp = chebpts(n); 
[t,r] = meshgrid(trigpts(m, [-pi, pi]), cp((n+1)/2:end));
F = f(t,r, 'polar');
G = sample(f, m, (n+1)/2);
pass(4) = norm(F(:) - G(:), inf) < tol;
[U, D, V] = sample(f, m, (n+1)/2);
G = U * D * V.';
pass(5) = norm(F(:) - G(:), inf) < tol;

% m odd and n odd
m = 31; 
n = 2*20+1;
cp = chebpts(n); 
[t,r] = meshgrid(trigpts(m, [-pi, pi]), cp((n+1)/2:end));
F = f(t,r, 'polar');
G = sample(f, m, (n+1)/2);
pass(6) = norm(F(:) - G(:), inf) < tol;
[U, D, V] = sample(f, m, (n+1)/2);
G = U * D * V.';
pass(7) = norm(F(:) - G(:), inf) < tol;

% m odd and n even
m = 31; 
n = 2*20-1;
cp = chebpts(n); 
[t,r] = meshgrid(trigpts(m, [-pi, pi]), cp((n+1)/2:end));
F = f(t,r, 'polar');
G = sample(f, m, (n+1)/2);
pass(8) = norm(F(:) - G(:), inf) < tol;
[U, D, V] = sample(f, m, (n+1)/2);
G = U * D * V.';
pass(9) = norm(F(:) - G(:), inf) < tol;

% m even and n even
m = 30; 
n = 2*20-1;
cp = chebpts(n); 
[t,r] = meshgrid(trigpts(m, [-pi, pi]), cp((n+1)/2:end));
F = f(t,r, 'polar');
G = sample(f, m, (n+1)/2);
pass(10) = norm(F(:) - G(:), inf) < tol;
[U, D, V] = sample(f, m, (n+1)/2);
G = U * D * V.';
pass(11) = norm(F(:) - G(:), inf) < tol;

% Sample should return all ones for the function 1.
f = diskfun(@(x,y) 1 + 0*x);
F = sample(f, 128, 128);
pass(12) = norm(F(:) - 1, inf) < tol;

% Check that errors are caught
try
    F = sample(f, 0, 20);
    pass(13) = false;
catch ME
    pass(13) = strcmp(ME.identifier, 'CHEBFUN:DISKFUN:sample:inputs');
end

try
    F = sample(f, 20, 0);
    pass(14) = false;
catch ME
    pass(14) = strcmp(ME.identifier, 'CHEBFUN:DISKFUN:sample:inputs');
end

end
function pass = test_constructor( ) 
% Test the spherefun constructor 

% Get tolerance:
tol = 2e3*chebfunpref().cheb2Prefs.chebfun2eps;

u = redefine_function_handle(@(x,y,z) x.^2 + y.^2 + z.^2);
v = u;
w = u;
f = spherefunv(u, v, w);
pass(1) = ( SampleError(f, u, v, w) < tol ); 

u = redefine_function_handle(@(x,y,z) exp(-cos(pi*(x+y+z))));
v = redefine_function_handle(@(x,y,z) x.*sin(x.*y));
w = redefine_function_handle(@(x,y,z) 1 + exp(x.*cos(x.*y)));
f = spherefunv(u, v, w);
pass(2) = ( SampleError(f, u, v, w) < tol );

u = redefine_function_handle(@(x,y,z) x.*y.*z);
uv = redefine_function_handle(@(x,y,z) x*y*z);
f = spherefunv(u, u, u);
g = spherefunv(uv, uv, uv, 'vectorize');
pass(3) = ( norm(f - g) < tol );

% Test errors
u = @(x,y,z) x;
try
    g = spherefunv(u);
    pass(4) = false;
catch ME
    pass(4) = strcmp(ME.identifier, 'CHEBFUN:SPHEREFUNV:spherefunv:arrayValued');
end

try
    g = spherefunv(u, u);
    pass(5) = false;
catch ME
    pass(5) = strcmp(ME.identifier, 'CHEBFUN:SPHEREFUNV:spherefunv:arrayValued');
end

try
    g = spherefunv(u, u, u, u);
    pass(6) = false;
catch ME
    pass(6) = strcmp(ME.identifier, 'CHEBFUN:SPHEREFUNV:spherefunv:arrayValued');
end

end

function f = redefine_function_handle(f)
% nargin(f) = 2, then we are already on the sphere, if nargin(f) = 3,
% then do change of variables:

if ( nargin(f) == 3 )
    % Wrap f so it can be evaluated in spherical coordinates
    f = @(lam, th) spherefun.sphf2cartf(f, lam, th, 0);
end

end

function sample_error = SampleError(h, u, v, w) 
m = 13; 
n = m;
[x, y] = getPoints(m, n);
[L2, T2] = meshgrid(x, y);
H = h(L2(:), T2(:));
U = [u(L2(:), T2(:)).'; 
    v(L2(:), T2(:)).'; 
    w(L2(:), T2(:)).'];
sample_error = norm(H(:) - U(:), inf);
end

function [x, y] = getPoints(m, n)

x = trigpts(2*n,[-pi pi]);
y = linspace(0, pi, m).';

end
function pass = test_constructor( ) 
% Test the diskfun constructor 

% Get tolerance:
tol = 2e3*chebfunpref().cheb2Prefs.chebfun2eps;

u = redefine_function_handle(@(x,y) x.^2 + y.^2);
v = u;
f = diskfunv(u, v, 'polar');
pass(1) = ( SampleError(f, u, v) < tol ); 
% try Cartesian
u2 = @(x,y) x.^2 + y.^2;
f = diskfunv(u2, u2);
pass(2) = ( SampleError(f, u, v) < tol );


u = redefine_function_handle(@(x,y) exp(-cos(pi*(x+y))));
v = redefine_function_handle(@(x,y) x.*sin(x.*y));
f = diskfunv(u, v, 'polar');
pass(3) = ( SampleError(f, u, v) < tol );

u = redefine_function_handle(@(x,y) x.*y);
uv = redefine_function_handle(@(x,y) x*y);
f = diskfunv(u, u);
g = diskfunv(uv, uv, 'vectorize');
pass(4) = ( norm(f - g) < tol );
%vectorize and polar 
u = @(t,r) r.*cos(t);
uv = @(t,r)  r*cos(t);
f = diskfunv(u, u, 'polar');
g = diskfunv(uv, uv, 'vectorize', 'polar');
pass(5) = ( norm(f - g) < tol );


% Test errors
u = @(x,y) x;
try
    g = diskfunv(u);
    pass(6) = false;
catch ME
    pass(6) = strcmp(ME.identifier, 'CHEBFUN:DISKFUNV:diskfunv:arrayValued');
end

try
    g = diskfunv(u, u, u);
    pass(7) = false;
catch ME
    pass(7) = strcmp(ME.identifier, 'CHEBFUN:DISKFUNV:diskfunv:arrayValued');
end

try 
g = diskfunv(u, 'polar');
    pass(8) = false;
catch ME
    pass(8) = strcmp(ME.identifier, 'CHEBFUN:DISKFUNV:diskfunv:arrayValued');
end
end


function sample_error = SampleError(h, u, v) 
m = 13; 
n = m;
[x, y] = getPoints(m, n);
[t, r] = meshgrid(x, y);
H = h(t(:), r(:), 'polar');
U = [u(t(:), r(:)).'; 
    v(t(:), r(:)).'];
sample_error = norm(H(:) - U(:), inf);
end

function [x, y] = getPoints(m, n)

x = trigpts(2*n, [-pi pi]);
y = chebpts(m);
y = y(ceil(m/2):end); 

end

function f = redefine_function_handle(f)
    % Wrap Cartesian f so it can be evaluated in polar coordinates
    
    f = @(th, r) diskfun.pol2cartf(f,th, r);
  
end
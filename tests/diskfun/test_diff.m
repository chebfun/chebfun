function pass = test_diff( )

tol = 2e2*chebfunpref().cheb2Prefs.chebfun2eps;

% Simple tests:
f = diskfun(@(th, r) r.^3.*(sin(th)+cos(th)), 'polar');  
fx = diff(f,1);
exact = @(th, r) r.^2.*(2*(cos(th)).^2 + 2*sin(th).*cos(th)+1);  
pass(1) = SampleError( exact, fx ) < tol;
fy = diff(f,2);
exact = @(th, r) r.^2.*(3*sin(th).^2+cos(th).^2+2*sin(th).*cos(th));  
pass(2) = SampleError( exact, fy ) < tol;


a = pi/4;  %constant derivatives
f = diskfun(@(th, r) r.*sin((th-a)), 'polar');  
fx = diff(f,1);
exact = @(th,r) -sin(a)+0*th;
pass(3) = SampleError( exact, fx ) < tol;
fy = diff(f,2);
exact = @(th,r) cos(a)+0*th;
pass(4) = SampleError( exact, fy ) < tol;

% Gaussian
g = @(x, y) 4/pi*exp(-4*((x-.3).^2+(y+.4).^2));
f = diskfun(g);
% dx using diff(f, 1)
fx = diff(f, 1);
exact = @(x,y) -32/pi*(x-.3).*exp(-4*((x-.3).^2+(y+.4).^2));
exact = redefine_function_handle(exact); 
pass(5) = SampleError( exact, fx ) < 2e1*tol;

% dx using diffx(f)
fx = diffx(f);
pass(6) = SampleError( exact, fx ) < 2e1*tol;

% dy using diff(f, 2)
fy = diff(f, 2);
exact = @(x,y) -32/pi*(y+.4).*exp(-4*((x-.3).^2+(y+.4).^2));
exact = redefine_function_handle(exact); 
pass(7) = SampleError( exact, fy ) < 2e1*tol;

% dy using diffy(f)
fy = diffy(f);
pass(8) = SampleError( exact, fy ) < 2e1*tol;

end
%%%%%%%%%%%%%%%%%%%%%

function sample_error = SampleError(h, g)
m = 7; 
n = m-1;
[x, y] = getPoints(m, n);
[L2, T2] = meshgrid(x, y);
F = h(L2, T2);
approx = fevalm(g, x, y);
sample_error = norm(F(:) - approx(:), inf);
end

function [x, y] = getPoints(m, n)

x = trigpts(2*n, [-pi pi]);
y = chebpts(m);
y = y(ceil(m/2):end); 

end

function f = redefine_function_handle(f)
    % Wrap Cartesian f so it can be evaluated in polar coordinates
    %if (coords=='cart')
    f = @(th, r) diskfun.pol2cartf(f,th, r);
    %end

end

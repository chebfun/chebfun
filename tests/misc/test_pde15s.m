function pass = test_pde15s(pref)

if ( nargin == 0 )
    pref = chebpref();
end

pass = 1;

return

%% pde15s_demos
% Some demos of pde15s (and pdesolve)

%% Advection
close all
d = [-1, 1]; x = chebfun('x', d);
u = exp(3*sin(pi*x));
plot(u); hold on
f = @(u, D) -D(u);
opts = pdeset('eps', 1e-4, 'abstol', 1e-4, 'reltol', 1e-4, 'plot', 1);
pde15s(f, 0:.05:2, u, 'periodic', opts);

%% Nonuniform Advection
close all
d = [-1, 1]; x = chebfun('x', d);
u = exp(3*sin(pi*x));
f = @(u, t, x, D) -(1+0.3*sin(pi*x)).*D(u) + 1e-6*D(u, 2);
opts = pdeset('eps', 1e-4, 'abstol', 1e-4, 'reltol', 1e-4, 'plot', 1);
pde15s(f, 0:.05:2, u, 'periodic', opts);

%% Advection-diffusion
close all
d = [-1, 1]; x = chebfun('x', d);
u = (1-x.^2).*exp(-30*(x+.5).^2);
f = @(u, t, x, D) -D(u)+.002*D(u, 2);
opts = pdeset('eps', 1e-4, 'abstol', 1e-4, 'reltol', 1e-4, 'plot', 1);
pde15s(f, 0:.05:2, u, 'dirichlet', opts);

%% Advection-diffusion2
close all
E = 1e-1;
d = [-3*pi/4, pi]; x = chebfun('x', d);
u = sin(2*x);
f = @(u, t, x, D) E*D(u, 2)+D(u);
lbc = {'neumann', 0};
rbc = {'dirichlet', 1};
bc = struct; bc.left = lbc; bc.right = rbc;
opts = pdeset('HoldPlot', 'on');
tt = linspace(0, 5, 21);
uu = pde15s(f, tt, u, bc, opts);

%% Advection-diffusion3 (Time depended rhs bc)
close all
E = 1e-1;
d = [-3*pi/4, pi]; x = chebfun('x', d);
u = sin(2*x);
f = @(u, t, x, D) E*D(u, 2)+D(u);
lbc = 'neumann';
rbc = @(u, t, x, D) u - .1*sin(t);
bc = struct; bc.left = lbc; bc.right = rbc;
opts = pdeset('holdPlot', 'on');
tt = linspace(0, 3, 101);
uu = pde15s(f, tt, u, bc, opts);

%% Allen-Cahn
close all
d = [-1, 1]; x = chebfun('x', d);
u = (1-x.^2).^2.*(1+sin(12*x))/2;
f = @(u, D) u.*(1-u) + 5e-4*D(u, 2);
pde15s(f, 0:0.1:5, u, 'neumann');

%% Allen-Cahn 
close all
d = [-1, 1]; x = chebfun('x', d);
u = .53*x-.47*sin(1.5*pi*x);
f = @(u, D) u.*(1-u.^2) + 5e-4*D(u, 2);
bc.left = struct('op', 'dirichlet', 'val', -1);
bc.right = struct('op', 'dirichlet', 'val', 1);
opts = pdeset('Ylim', [-1.1 1.1]);
pde15s(f, 0:0.1:5, u, bc, opts);

%% Burgers
close all
d = [-1, 1]; x = chebfun('x', d);
u = (1-x.^2).*exp(-30*(x+.5).^2);
f = @(u, D) -D(u.^2)+.01*D(u, 2);
pde15s(f, 0:.05:3, u, 'dirichlet');

%% KS
close all
d = [-1, 1]; x = chebfun('x', d);
I = eye(d); D = diff(d);
u = 1 + 0.5*exp(-40*x.^2);
bc = struct;
bc.left = @(u, D) [u-1 ; D(u)];
bc.right = @(u, D) [u-1 ; D(u)];
f = @(u, D) u.*D(u)-D(u, 2)-0.006*D(u, 4);
u = pde15s(f, 0:.01:.5, u, bc);

%% Cahn-Hilliard - not working!
close all
E = 1e-1;
d = [-1, 1]; x = chebfun('x', d);
u = cos(pi*x)-exp(-6*pi*x.^2);
plot(u)
opts = pdeset('eps', 1e-4, 'abstol', 1e-4, 'reltol', 1e-4, 'plot', 1);
bc = struct;
bc.left = @(u, D) [u+1 ; D(u)];
bc.right = @(u, D) [u+1 ; D(u)];
f = @(u, t, x, D) -D(u, 4) + D(u.^3, 2)-D(u, 2);
tt = linspace(0, .001, 101);
uu = pde15s(f, tt, u, bc, opts);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% A selection of PDE system examples
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%

close all
d = [-1, 1]; x = chebfun('x', d);
u = [ chebfun(1, d)  chebfun(1, d) ];
f = @(u, v, t, x, D) [ -u + (x + 1).*v + 0.1*D(u, 2) , ...
                        u - (x + 1).*v + 0.2*D(v, 2) ];
bc.left = @(u, v, D) [D(u), D(v)];  bc.right = bc.left;   % New way
opts = pdeset('plot', 1, 'N', 32);
uu = pde15s(f, 0:.05:3, u, bc, opts);


%%

close all
d = [-1, 1]; x = chebfun('x', d);
u = [ 1-erf(10*(x+0.7)) , 1 + erf(10*(x-0.7)) , chebfun(0, d) ];
f = @(u, v, w, D)      [ 0.1*D(u, 2) - 100*u.*v , ...
                      0.2*D(v, 2) - 100*u.*v , ...
                     .001*D(w, 2) + 2*100*u.*v ];                  
bc = 'neumann';     
uu = pde15s(f, 0:.1:3, u, bc);

%%

close all
% Crazy nonlinear boundary conditions
d = [-1, 1]; x = chebfun('x', d);
u = [ chebfun(1, d)  chebfun(1, d) ];
f = @(u, v, t, x, D) [ -u + (x + 1).*v + 0.1*D(u, 2) , ...
                    u - (x + 1).*v + 0.2*D(v, 2) ];
bc.left =  @(u, v, t, x, D) [D(u)+t*sin(u)./v,  D(v)];
bc.right = @(u, v, t, x, D) [D(u),              D(v).*v+sin(5*pi*t)];
uu = pde15s(f, 0:.05:1, u, bc);

%%

close all
% Maxwell's Equations
d = [-1, 1]; x = chebfun('x', d);
u = exp(-20*x.^2) .* sin(14*x);  u = [u -u];
D = diff(d); I = eye(d); Z = zeros(d);
f = @(u, v, D) [D(v) ; D(u)];
bc.left = @(u, v, D) u; bc.right = @(u, v, D) v;        % New way
opt = pdeset('eps', 1e-6, 'Ylim', pi/2*[-1 1], 'AbsTol', 1e-6, 'RelTol', 1e-6);
uu = pde15s(f, 0:.05:2, u, bc, opt, 64);

end

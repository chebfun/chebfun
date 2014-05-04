function pass = test_pde15s(pref, flag)

if ( nargin == 0 )
    pref = chebfunpref();
end

pass = 1;

if ( nargin < 2  )
    pause(.1+ .1*rand())
else
    try
        doTests()
    catch
        pass = 0;
    end
end

end

function doTests()

%% pde15s_demos
% Some demos of pde15s (and pdesolve)

%% Advection
close all
d = [-1, 1]; x = chebfun('x', d);
u = exp(3*sin(pi*x));
plot(u); hold on
f = @(u) -diff(u);
opts = pdeset('eps', 1e-4, 'abstol', 1e-4, 'reltol', 1e-4, 'plot', 1);
uu = pde15s(f, 0:.05:2, u, 'periodic', opts);

%% Nonuniform Advection
close all
d = [-1, 1]; x = chebfun('x', d);
u = exp(3*sin(pi*x));
f = @(t, x, u) -(1+0.3*sin(pi*x)).*diff(u) + 1e-6*diff(u, 2);
opts = pdeset('eps', 1e-4, 'abstol', 1e-4, 'reltol', 1e-4, 'plot', 1);
uu = pde15s(f, 0:.05:2, u, 'periodic', opts);

%% Advection-diffusion
close all
d = [-1, 1]; x = chebfun('x', d);
u = (1-x.^2).*exp(-30*(x+.5).^2);
f = @(t, x, u) -diff(u)+.002*diff(u, 2);
opts = pdeset('eps', 1e-4, 'abstol', 1e-4, 'reltol', 1e-4, 'plot', 1);
uu = pde15s(f, 0:.05:2, u, 'dirichlet', opts);
waterfall(uu), shg

%% Advection-diffusion2
close all
E = 1e-1;
d = [-3*pi/4, pi]; x = chebfun('x', d);
u = sin(2*x);
f = @(t, x, u) E*diff(u, 2)+diff(u);
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
f = @(t, x, u) E*diff(u, 2)+diff(u);
lbc = 'neumann';
rbc = @(t, x, u) u - .1*sin(t);
bc = struct; bc.left = lbc; bc.right = rbc;
opts = pdeset('holdPlot', 'on');
tt = linspace(0, 3, 51);
uu = pde15s(f, tt, u, bc, opts);

%% Allen-Cahn
close all
d = [-1, 1]; x = chebfun('x', d);
u = (1-x.^2).^2.*(1+sin(12*x))/2;
f = @(u) u.*(1-u) + 5e-4*diff(u, 2);
pde15s(f, 0:0.1:5, u, 'neumann');

%% Allen-Cahn 
close all
d = [-1, 1]; x = chebfun('x', d);
u = .53*x-.47*sin(1.5*pi*x);
f = @(u) u.*(1-u.^2) + 5e-4*diff(u, 2);
bc.left = struct('op', 'dirichlet', 'val', -1);
bc.right = struct('op', 'dirichlet', 'val', 1);
opts = pdeset('Ylim', [-1.1 1.1]);
pde15s(f, 0:0.1:5, u, bc, opts);

%% Burgers
close all
d = [-1, 1]; x = chebfun('x', d);
u = (1-x.^2).*exp(-30*(x+.5).^2);
f = @(u) -diff(u.^2)+.01*diff(u, 2);
pde15s(f, 0:.05:3, u, 'dirichlet');

%% KS
close all
d = [-1, 1]; x = chebfun('x', d);
u = 1 + 0.5*exp(-40*x.^2);
bc = struct;
bc.left = @(u) [u-1 ; diff(u)];
bc.right = @(u) [u-1 ; diff(u)];
f = @(u) u.*diff(u)-diff(u, 2)-0.006*diff(u, 4);
u = pde15s(f, 0:.01:.5, u, bc);

%% Cahn-Hilliard - not working!
close all
E = 1e-1;
d = [-1, 1]; x = chebfun('x', d);
u = cos(pi*x)-exp(-6*pi*x.^2);
plot(u)
opts = pdeset('eps', 1e-4, 'abstol', 1e-4, 'reltol', 1e-4, 'plot', 1);
bc = struct;
bc.left = @(u) [u+1 ; diff(u)];
bc.right = @(u) [u+1 ; diff(u)];
f = @(t, x, u) -diff(u, 4) + diff(u.^3, 2)-diff(u, 2);
tt = linspace(0, .001, 101);
uu = pde15s(f, tt, u, bc, opts);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% A selection of PDE system examples
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%

close all
d = [-1, 1]; x = chebfun('x', d);
u = [ chebfun(1, d)  chebfun(1, d) ];
f = @(t, x, u, v) [ -u + (x + 1).*v + 0.1*diff(u, 2) , ...
                     u - (x + 1).*v + 0.2*diff(v, 2) ];
bc.left = @(u, v) [diff(u), diff(v)];  bc.right = bc.left;   % New way
opts = pdeset('plot', 1, 'N', 32);
uu = pde15s(f, 0:.05:3, u, bc, opts);


%%

close all
d = [-1, 1]; x = chebfun('x', d);
u = [ 1-erf(10*(x+0.7)) , 1 + erf(10*(x-0.7)) , chebfun(0, d) ];
f = @(u, v, w)      [ 0.1*diff(u, 2) - 100*u.*v , ...
                      0.2*diff(v, 2) - 100*u.*v , ...
                     .001*diff(w, 2) + 2*100*u.*v ];                  
bc = 'neumann';     
uu = pde15s(f, 0:.1:3, u, bc);

%%

close all
% Crazy nonlinear boundary conditions
d = [-1, 1]; x = chebfun('x', d);
u = [ chebfun(1, d)  chebfun(1, d) ];
f = @(t, x, u, v) [ -u + (x + 1).*v + 0.1*diff(u, 2) , ...
                     u - (x + 1).*v + 0.2*diff(v, 2) ];
bc = struct;
bc.left =  @(t, x, u, v) [diff(u)+t*sin(u)./v,  diff(v)];
bc.right = @(t, x, u, v) [diff(u),              diff(v).*v+sin(5*pi*t)];
uu = pde15s(f, 0:.05:1, u, bc);

%%

close all
% Maxwell's Equations
d = [-1, 1]; x = chebfun('x', d);
u = exp(-20*x.^2) .* sin(14*x);  u = [u -u];
f = @(u, v) [diff(v) ; diff(u)];
bc.left = @(u, v) u; bc.right = @(u, v) v;        % New way
opt = pdeset('eps', 1e-6, 'Ylim', pi/2*[-1 1], 'AbsTol', 1e-6, 'RelTol', 1e-6);
uu = pde15s(f, 0:.05:2, u, bc, opt, 64);

%% 
%chebop-style synatx

close all
x = chebfun('x');
u = 1 + 0.5*exp(-40*x.^2);
bcc = @(x, u) [u(-1)-1 ; feval(diff(u),-1) ; u(1)-1 ; feval(diff(u),1)];
f = @(u) u.*diff(u) - diff(u, 2) - 0.006*diff(u, 4);
opts = pdeset('Ylim', [-30 30], 'PlotStyle', {'LineWidth', 2});
uu = pde15s(f, 0:.025:.5, u, bcc, opts);
surf(uu, 0:.025:.5)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% V4 syntax
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Example 1: Nonuniform advection
  x = chebfun('x',[-1 1]);
  u = exp(3*sin(pi*x));
  f = @(u,t,x,diff) -(1+0.6*sin(pi*x)).*diff(u);
  uu = pde15s(f,0:.05:.5,u,'periodic');

%% Example 2: Kuramoto-Sivashinsky
  d = domain(-1,1);
  x = chebfun('x');
  I = eye(d); D = diff(d);
  u = 1 + 0.5*exp(-40*x.^2);
  bc.left = struct('op',{I,D},'val',{1,2});
  bc.right = struct('op',{I,D},'val',{1,2});
  f = @(u,diff) u.*diff(u)-diff(u,2)-0.006*diff(u,4);
  try
    uu = pde15s(f,0:.01:.5,u,bc);
  catch ME
      if ( ~strcmp(ME.identifier, 'CHEBFUN:pde15s:bcstruct') )
          rethrow(ME)
      end
  end

%% Example 3: Chemical reaction (system)
   x = chebfun('x',[-1 1]);  
   u = [ 1-erf(10*(x+0.7)) , 1 + erf(10*(x-0.7)) , 0 ];
   f = @(u,v,w,diff)  [ .1*diff(u,2) - 100*u.*v , ...
                        .2*diff(v,2) - 100*u.*v , ...
                        .001*diff(w,2) + 2*100*u.*v ];
   bc = 'neumann';     
   uu = pde15s(f,0:.1:3,u,bc);

end

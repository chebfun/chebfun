function pass = test_pde15s(pref, varargin)

if ( nargin == 0 )
    pref = chebfunpref();
end
tol = 1e5*pref.chebfuneps;


%% Test chebtech1 works (see #910)
% Get global prefs:
globalPref = chebfunpref(); 
try

    % Set 1st kind grids:
    chebfunpref.setDefaults('tech', @chebtech1);
    f = chebfun(@(x) sin(pi*x), [-2.5 3]);
    bc.left = @(u) diff(u);
    bc.right = 0;
    opts = pdeset('plot','off');
    uu = pde15s(@(t,x,u) .1*diff(u,2) + diff(u), 0:.1:1, f, bc, opts);

    % Set 2nd kind grids:
    chebfunpref.setDefaults('tech', @chebtech2);
    f = chebfun(@(x) sin(pi*x), [-2.5 3]);
    bc.left = @(u) diff(u);
    bc.right = 0;
    opts = pdeset('plot','off');
    vv = pde15s(@(t,x,u) .1*diff(u,2) + diff(u), 0:.1:1, f, bc, opts);

    % Go Compare:
    pass(1) = norm( uu - vv ) < tol ; 
    
    % Reset preferences: 
    chebfunpref.setDefaults(globalPref); 
    
catch ME
    % Reset preferences: 
    chebfunpref.setDefaults(globalPref); 
    rethrow(ME)
end

if ( nargin == 2)
    doTests();
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
uu = pde23t(f, 0:.05:2, u, 'periodic', opts);

%% Nonuniform Advection
close all
d = [-1, 1]; x = chebfun('x', d);
u = exp(3*sin(pi*x));
f = @(t, x, u) -(1+0.3*sin(pi*x)).*diff(u);
opts = pdeset('eps', 1e-4, 'abstol', 1e-4, 'reltol', 1e-4, 'plot', 1);
uu = pde23t(f, 0:.05:2.5, u, 'periodic', opts);

%% Advection-diffusion
close all
d = [-1, 1]; x = chebfun('x', d);
u = (1-x.^2).*exp(-30*(x+.5).^2);
f = @(t, x, u) -diff(u)+.002*diff(u, 2);
opts = pdeset('plot', 1);
uu = pde15s(f, 0:.05:2, u, 'dirichlet', opts);
waterfall(uu), shg

%% Advection-diffusion2
close all
E = 1e-1;
d = [-3*pi/4, pi]; x = chebfun('x', d);
u = sin(2*x);
f = @(t, x, u) E*diff(u, 2)+diff(u);
lbc = {'neumann', 0};
rbc = {'dirichlet', 0};
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
rbc = @(t, u) u - .1*sin(t);
bc = struct; bc.left = lbc; bc.right = rbc;
opts = pdeset('holdPlot', 'on');
tt = linspace(0, 3, 51);
uu = pde15s(f, tt, u, bc, opts);

%% Advection-diffusion4 (Time depended rhs bc)
close all
E = 1e-1;
d = [-3*pi/4, pi]; x = chebfun('x', d);
u = sin(2*x);
f = @(t, x, u) E*diff(u, 2)+diff(u);
rbc = @(t, u) u - .1*sin(t);
mid = @(t, x, u) feval(u, d(1)) - 1;
bc = struct; bc.left = []; bc.right = rbc; bc.middle = mid;
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
bc = struct;
% bc.left = struct('op', 'dirichlet', 'val', -1);
% bc.right = struct('op', 'dirichlet', 'val', 1);
bc.left = @(u) u + 1;
bc.right = @(u) u - 1;
opts = pdeset('Ylim', [-1.1 1.1]);
pde15s(f, 0:0.1:3, u, bc, opts);

%% Burgers
close all
d = [-1, 1]; x = chebfun('x', d);
u = (1-x.^2).*exp(-30*(x+.5).^2);
f = @(u) -diff(u.^2)+.01*diff(u, 2);
pde15s(f, 0:.1:3, u, 'dirichlet');

%% KS
close all
d = [-1, 1]; x = chebfun('x', d);
u = 1 + 0.5*exp(-40*x.^2);
bc = struct;
bc.left = @(u) [u-1 ; diff(u)];
bc.right = @(u) [u-1 ; diff(u)];
f = @(u) u.*diff(u)-diff(u, 2)-0.006*diff(u, 4);
u = pde15s(f, 0:.01:.5, u, bc);

%% Cahn-Hilliard
close all
E = 1e-1;
d = [-1, 1]; x = chebfun('x', d);
u = cos(pi*x)-exp(-6*pi*x.^2);
plot(u)
opts = pdeset('plot', 1);
bc = struct;
bc.left = @(u) [u+1 ; diff(u)];
bc.right = @(u) [u+1 ; diff(u)];
f = @(t, x, u) -diff(u, 4) + diff(u.^3, 2)-diff(u, 2);
tt = linspace(0, .02, 51);
uu = pde15s(f, tt, u, bc, opts);

%% integral operator:

dom = [-1 1];
t = 0:.1:4;
pdefun = @(t,x,u) .02.*diff(u,2)+cumsum(u).*sum(u);
bc.left = 'dirichlet';
bc.right = 'dirichlet';
x = chebfun(@(x) x, dom);
u0 = (1-x.^2).*exp(-30.*(x+.5).^2);
opts = pdeset('Eps',1e-6,'Ylim',[0,1.4]);
[t, u] = pde15s(pdefun, t, u0, bc, opts);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% A selection of PDE system examples
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%

close all
d = [-1, 1]; x = chebfun('x', d);
u = [ chebfun(1, d) ; chebfun(1, d) ];
f = @(t, x, u, v) [ -u + (x + 1).*v + 0.1*diff(u, 2) ; ...
                     u - (x + 1).*v + 0.2*diff(v, 2) ];
bc.left = @(u, v) [diff(u), diff(v)];  bc.right = bc.left;   % New way
opts = pdeset('plot', 1);
uu = pde15s(f, 0:.05:3, u, bc, opts);


%%

close all
d = [-1, 1]; x = chebfun('x', d);
u = [ 1-erf(10*(x+0.7)) ; 1 + erf(10*(x-0.7)) ; chebfun(0, d) ];
f = @(u, v, w)      [ 0.1*diff(u, 2) - 100*u.*v ; ...
                      0.2*diff(v, 2) - 100*u.*v ; ...
                     .001*diff(w, 2) + 2*100*u.*v ];                  
bc = 'neumann';     
uu = pde15s(f, 0:.1:3, u, bc);

%%

close all
% Crazy nonlinear boundary conditions
d = [-1, 1]; x = chebfun('x', d);
u = [ chebfun(1, d) ; chebfun(1, d) ];
f = @(t, x, u, v) [ -u + (x + 1).*v + 0.1*diff(u, 2) ; ...
                     u - (x + 1).*v + 0.2*diff(v, 2) ];
bc = struct;
bc.left =  @(t, x, u, v) [diff(u)+t*sin(u)./v ;  diff(v)];
bc.right = @(t, x, u, v) [diff(u) ;              diff(v).*v+sin(5*pi*t)];
uu = pde15s(f, 0:.05:1, u, bc);

%%

close all
% Maxwell's Equations
d = [-1, 1]; x = chebfun('x', d);
u = exp(-20*x.^2) .* sin(14*x);  u = [u ; -u];
f = @(u, v) [diff(v) ; diff(u)];
bc.left = @(u, v) u; bc.right = @(u, v) v;        % New way
opt = pdeset('eps', 1e-6, 'Ylim', pi/2*[-1 1], 'AbsTol', 1e-6, 'RelTol', 1e-6);
[t, u, v] = pde23t(f, 0:.05:2, u, bc, opt);

%% 
%chebop-style synatx

close all
x = chebfun('x');
u = 1 + 0.5*exp(-40*x.^2);
bcc = @(x,t, u) [u(-1)-1 ; feval(diff(u),-1) ; u(1)-1 ; feval(diff(u),1)];
f = @(u) u.*diff(u) - diff(u, 2) - 0.006*diff(u, 4);
opts = pdeset('Ylim', [-30 30], 'PlotStyle', {'LineWidth', 2});
uu = pde15s(f, 0:.025:.5, u, bcc, opts);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% V4 syntax
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Example 1: Nonuniform advection
    close all
    d = [-1, 1]; x = chebfun('x', d);
    u = exp(3*sin(pi*x));
    f = @(u,t,x,diff) -(1+0.3*sin(pi*x)).*diff(u);
    opts = pdeset('eps', 1e-4, 'abstol', 1e-4, 'reltol', 1e-4, 'plot', 1);
    uu = pde23t(f, 0:.05:2, u, 'periodic', opts);

%% Example 2: Kuramoto-Sivashinsky
  close all
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
      if ( ~strcmp(ME.identifier, 'CHEBFUN:CHEBFUN:pde15s:bcstruct') )
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
   [t, u, v, w] = pde15s(f,0:.1:3,u,bc);
   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MISC
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
%% A periodic example which was failing (#2349)
    dom = [0, 1];
    t = [0, .1];
    bc = 'periodic';
    pdefun = @(t,x,u,v) [diff(v); diff(u)];
    x = chebfun(@(x) x, dom);
    ic = [sin(2*pi*x), 0];
    [t, u, v] = pde23t(pdefun, t, ic, bc);   
   
end

function pass = test_minres(pref) 
% Test preconditioned minres.

% Tolerance: 
if ( nargin == 0 )
    pref = cheboppref;
end
tol = 1e2*pref.bvpTol;

%% Example 1: (-u_xx=f, bc=0, sum(f)==0)
a = @(x) 1+0*x;
c = @(x) -0*x;
L = @(x,u) -diff(a(x).*diff(u)) + c(x).*u;

% Chebop solve:
N = chebop( L );
N.bc = 0; 
f = chebfun(@(x) 1-3*x.^2);
u = N \ f; 

v = minres(N, f);

% Errors:
pass(1) = ( norm( u - v ) < tol );

%% Example 2: (-(a(x)*u_x)_x+c(x)*u=f, bc=0, sum(f)==0)
a = @(x) 2+cos(pi*x);
c = @(x) -10*x;
L = @(x,u) -diff(a(x).*diff(u)) + c(x).*u;

% Chebop solve:
N = chebop( L );
N.bc = 0; 
f = chebfun(@(x) 1-3*x.^2);
u = N \ f; 

v = minres(N, f);

% Errors:
pass(2) = ( norm( u - v ) < tol );

%% Example 3: (-(a(x)*u_x)_x + c(x)*u=f, bc=0, sum(f)==0)

a = @(x) 2+cos(pi*x);
c = @(x) 1-10*x.^2;
L = @(x,u) -diff(a(x).*diff(u)) + c(x).*u;

% Chebop solve:
N = chebop( L );
N.bc = 0; 
f = chebfun(@(x) 1-3*x.^2);
u = N \ f; 

v = minres(N, f);

% Errors:
pass(3) = ( norm( u - v ) < tol );

%% Example 4: (-(a(x)*u_x)_x + c(x)*u=f, bc=0, mean(f)~=0)
a = @(x) 2+cos(pi*x);
c = @(x) 1+10*x.^2;
L = @(x,u) -diff(a(x).*diff(u)) + c(x).*u;

% Chebop solve:
N = chebop( L );
N.bc = 0; 
f = chebfun(@(x) 1-2*x.^2);
u = N \ f; 

% Chebfun-ready pcg:
v = minres(N, f);

pass(4) = ( norm(u - v) < tol );

%% Example 5: (-(a(x)*u_x)_x + c(x)*u=f, bc~=0, mean(f)~=0)
a = @(x) 2+cos(pi*x);
c = @(x) 1-100*x.^2;
L = @(x,u) -diff(a(x).*diff(u)) + c(x).*u;

% Chebop solve:
N = chebop( L );
N.lbc = 1; 
N.rbc = -1; 
f = chebfun(@(x) 1-2*x.^2);
u = N \ f; 

% Chebfun-ready pcg:
v = minres(N, f, [], 40);

pass(5) = ( norm(u - v) < tol );

if ( all(pass) == 1 )
    pass = 1;
end

end
clear classes
clc

% We solve y'' + y' + x*y = 1, y(-1)=1, y(1)=0, using integral formulation.
% Let v = y'', so that y'=C*v + a, y = C*C*v + a*(x+1) + b, where C is the
% cumsum operator. We have
%
%    v + C*v + a + x*(C*C*v) + x*(a*(x+1)) + x*b = 1
%
%    y = 1 at left
%
%    y = 0 at right
%
% We solve as a problem for v, a, and b.

%%
% Basic operators.
I = linop.eye;
C = linop.cumsum;
Z = linop.zeros;

x = chebfun('x');
one = chebfun(1);
X = linop.diag(x);

%%
% This opeartor links v to y. It'll be handy a couple of times.
v2y = [C^2 x+1 one];

%%
% This operator applies to the chebmatrix {v;a;b}.
L = [ I + C + X*C^2, 1+x.*(x+1), x ];

%%
% Add on the boundary conditions and define the RHS. 
L = lbc(L,v2y,1);
L = rbc(L,v2y,0);
f = chebmatrix({one});

%%
% Here's what the discrete linear system looks like at n=6.
[A,b] = linSystem(L,f,6,@colloc2)

%%
% Solve for the second derivative.
v = L\f
plot(v{1}), ylabel('y''''')

%%
% Convert to the function y.
y = v2y*v;
plot(y{1}), ylabel('y')

%%
% Check residual.
y = y{1};
ShouldBeZero = norm( diff(y,2)+diff(y)+x.*y-one )

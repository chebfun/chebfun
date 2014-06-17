function pass = test_ode45( pref ) 
% A large battery of functions. 

if ( nargin < 1 ) 
    pref = chebfunpref; 
end 
tol = 100 * pref.eps; 

% Trivial check for ODE45 for now:
t0 = 0; Tend = 50;
F = chebfun2v(@(h,hp)hp, @(h,hp)-1-.01*hp, @(h,hp)1+0*h,[0 30 0 2]);
u0 = [2 0 0];  % u0 = [h(0) h'(0) x(0)]
options = odeset('RelTol',100*eps, 'events', @p1Event1);
[T, Y] = ode45(F, [t0, Tend], u0, options);
pass(1) = ( norm( Y(0,1) - u0(1) ) < tol );
end

function [value, isTerminal, direction] = p1Event1(t,y)
% Event handling function. Stop Chebfun2's ODE45 when the ball bounces.
value = y(1) - y(3);
isTerminal = 1;
direction = 0;
end
function pass = test_ode45( pref ) 
% A large battery of functions. 

if ( nargin < 1 ) 
    pref = chebfunpref; 
end 
tol = 100 * pref.eps;                % TODO: This isn't used.

% Trivial check for ODE45 for now:
t0 = 0; Tend = 30;
F = chebfun2v(@(h,hp)hp, @(h,hp)-1-.01*hp, @(h,hp)1+0*h,[0 30 0 2]);
u0 = [2 0 0];  % u0 = [h(0) h'(0) x(0)]
options = odeset('RelTol',100*eps, 'events', @p1Event1);
[T, Y] = ode45(F, [t0, Tend], u0, options);
pass(1) = ( norm( Y(0,1) - u0(1) ) < 1e-3 ); % TODO: This weems loose.

% The following fails if ode45 does not ensure the trajectory staying the
% domain of the chebfun2v. 
pass(2) = 1; 
try 
g = chebfun2v(@(x,y) x, @(x,y) y);   % identity chebfun2v
A = [2 -2 ; 0 1];                    % matrix of the system and eigenvalues
G = A*g; T = [0 1];                  % phase plane and time interval

initvals = [.1 .05; -.1 -.05; -.1,-.05; -.1,0; .1,0];
for k = 1:size(initvals,1)
    [ignored, y] = ode45(G,T,initvals(k,:));
end
initvals = [.1 .1; -.1 -.1];
for k = 1:size(initvals,1)
    [ignored, y] = ode45(G,2*T/3,initvals(k,:));
end
catch 
   pass(2) = 0;  
end

end

function [value, isTerminal, direction] = p1Event1(t, y)
% Event handling function. Stop Chebfun2's ODE45 when the ball bounces.
value = y(1) - y(3);
isTerminal = 1;
direction = 0;
end
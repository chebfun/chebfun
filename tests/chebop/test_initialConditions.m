
function pass = test_initialConditions(pref)
%TEST_INITIALCONDITIONS   Checks conditions for IVP/FVP problems.

%% Setup
if ( nargin == 0 )
    pref = cheboppref();
end
dom = [0, 1];
% [TODO]: This used to be 1e-12 until we had to loosen the absTol of the ODE
% solvers. Once/if that gets resolved, can we tighten the tolerance again?
tol = 5e-9;

%% First order IVP
N = chebop(@(x,u) diff(u) + sin(u), dom);
N.lbc = @(u) u - .5;
u = N\0;
pass(1) = abs(u(dom(1)) - .5) < tol;

%% First order FVP
N = chebop(@(x,u) diff(u) + sin(u), dom);
N.rbc = @(u) u - .5;
u = N\0;
pass(2) = abs(u(dom(end)) - .5) < tol;

%% Second order IVP
N = chebop(@(x,u) diff(u, 2) + sin(u), dom);
N.lbc = @(u) [u - .5; diff(u) - 1.3];
u = N\0;
up = diff(u);
pass(3) = abs(u(dom(1)) - .5) + abs(up(dom(1)) - 1.3) < tol;

%% Second order FVP
N = chebop(@(x,u) diff(u, 2) + sin(u), dom);
N.rbc = @(u) [u - .5; diff(u) - 1.3];
u = N\0;
up = diff(u);
pass(4) = abs(u(dom(end)) - .5) + abs(up(dom(end)) - 1.3) < tol;

%% First order IVP, no x
N = chebop(@(u) diff(u) + sin(u), dom);
N.lbc = @(u) u - .5;
u = N\0;
pass(5) = abs(u(dom(1)) - .5) < tol;

%% First order FVP, no x
N = chebop(@(u) diff(u) + sin(u), dom);
N.rbc = @(u) u - .5;
u = N\0;
pass(6) = abs(u(dom(end)) - .5) < tol;

%% Second order IVP, no x
N = chebop(@(u) diff(u, 2) + sin(u), dom);
N.lbc = @(u) [u - .5; diff(u) - 1.3];
u = N\0;
up = diff(u);
pass(7) = abs(u(dom(1)) - .5) + abs(up(dom(1)) - 1.3) < tol;

%% Second order FVP, no x
N = chebop(@(u) diff(u, 2) + sin(u), dom);
N.rbc = @(u) [u - .5; diff(u) - 1.3];
u = N\0;
up = diff(u);
pass(8) = abs(u(dom(end)) - .5) + abs(up(dom(end)) - 1.3) < tol;

%% First order coupled IVP
N = chebop(@(x,u, v, w) [diff(u) + v; diff(v) + w; diff(w) + u], dom);
N.lbc = @(u,v,w) [u - .5; v - 1.5; w + 2.5];
uvw = N\[0; 1; 2];
u = uvw{1}; v = uvw{2}; w = uvw{3};
pass(9) = abs(u(dom(1)) - .5) + abs(v(dom(1)) - 1.5) + abs(w(dom(1)) + 2.5) < tol;

%% First order coupled FVP
N = chebop(@(x,u, v, w) [diff(u) + v; diff(v) + w; diff(w) + u], dom);
N.rbc = @(u,v,w) [u - .5; v - 1.5; w + 2.5];
uvw = N\[0; 1; 2];
u = uvw{1}; v = uvw{2}; w = uvw{3};
pass(10) = abs(u(dom(end)) - .5) + abs(v(dom(end)) - 1.5) + abs(w(dom(end)) + 2.5) < tol;

%% Third order coupled IVP
N = chebop(@(x,u, v, w) [diff(u, 3) + v; diff(v, 3) + w; diff(w, 3) + u], dom);
N.lbc = @(u,v,w) [u - .5; v - 1.5; w + 2.5; ...
    diff(u) - 1; diff(v) - 2; diff(w) + 1; ...
    diff(u,2); diff(v,2) - .5; diff(w,2) + .5];
uvw = N\[0; 1; 2];
u = uvw{1}; v = uvw{2}; w = uvw{3};
up = diff(u); vp = diff(v); wp = diff(w);
upp = diff(u, 2); vpp = diff(v, 2); wpp = diff(w, 2);
pass(11) = abs(u(dom(1)) - .5) + abs(v(dom(1)) - 1.5) + abs(w(dom(1)) + 2.5) + ...
    abs(up(dom(1)) - 1) + abs(vp(dom(1)) - 2) + abs(wp(dom(1)) + 1) + ...
    abs(upp(dom(1))) + abs(vpp(dom(1)) - .5) + abs(wpp(dom(1)) + .5) < 100*tol;

%% Third order coupled FVP
N = chebop(@(x,u, v, w) [diff(u, 3) + v; diff(v, 3) + w; diff(w, 3) + u], dom);
N.rbc = @(u,v,w) [u - .5; v - 1.5; w + 2.5; ...
    diff(u) - 1; diff(v) - 2; diff(w) + 1; ...
    diff(u,2); diff(v,2) - .5; diff(w,2) + .5];
uvw = N\[0; 1; 2];
u = uvw{1}; v = uvw{2}; w = uvw{3};
up = diff(u); vp = diff(v); wp = diff(w);
upp = diff(u, 2); vpp = diff(v, 2); wpp = diff(w, 2);
pass(12) = abs(u(dom(end)) - .5) + abs(v(dom(end)) - 1.5) + abs(w(dom(end)) + 2.5) + ...
    abs(up(dom(end)) - 1) + abs(vp(dom(end)) - 2) + abs(wp(dom(end)) + 1) + ...
    abs(upp(dom(end))) + abs(vpp(dom(end)) - .5) + abs(wpp(dom(end)) + .5) < 15*tol;

%% Incorrect number of conditions

% Too few conditions.
N = chebop(@(x,u) diff(u, 2) + sin(u), dom);
N.lbc = @(u) u;
try
    mldivide(N, 0);
catch ME
    pass(13) = strcmp(ME.identifier, 'CHEBFUN:CHEBOP:solveivp:numConditions');
end

% Too many conditions
N = chebop(@(x,u) diff(u) + sin(u), dom);
N.lbc = @(u) [u; diff(u)];
try
    mldivide(N, 0);
catch ME
    pass(14) = strcmp(ME.identifier, 'CHEBFUN:CHEBOP:solveivp:numConditions');
end

end
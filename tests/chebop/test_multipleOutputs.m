function pass = test_multipleOutputs(~)
% TEST_MULTIPLEOUTPUTS   Check that syntax [u,v,w] = N\0 works.

%% BVP scalar case
dom  = [0, 2];
N = chebop(@(x,u) diff(u,2) + 0.01*sin(u), dom);
N.lbc = 0;
N.rbc = 2;
u = N\0;
pass(1) = isa(u, 'chebfun');

%% BVP scalar case with info
dom  = [0, 2];
N = chebop(@(x,u) diff(u,2) + 0.01*sin(u), dom);
N.lbc = 0;
N.rbc = 2;
[u, info] = N\0;
pass(2) = isa(u, 'chebfun') && isstruct(info);

%% IVP scalar case
dom  = [0, 2];
N = chebop(@(x,u) diff(u,2) + sin(u), dom);
N.lbc = @(u) [u; diff(u)-1];
u = N\0;
pass(3) = isa(u, 'chebfun');

%% IVP scalar case with info
dom  = [0, 2];
N = chebop(@(x,u) diff(u,2) + sin(u), dom);
N.lbc = @(u) [u; diff(u)-1];
[u, info] = N\0;
pass(4) = isa(u, 'chebfun') && isstruct(info);

%% BVP two var case, one output
dom  = [0, 2];
N = chebop(@(x,u,v) [diff(u,2) + 0.1*sin(u) + v; diff(v,2) + 0.1*cos(v) + u], dom);
N.lbc = @(u,v) [u - 1 ; diff(v)];
N.rbc = @(u,v) [diff(v); v + 1];
u = N\0;
pass(5) = isa(u, 'chebmatrix');

%% BVP two var case, two outputs
dom  = [0, 2];
N = chebop(@(x,u,v) [diff(u,2) + 0.1*sin(u) + v; diff(v,2) + 0.1*cos(v) + u], dom);
N.lbc = @(u,v) [u - 1 ; diff(v)];
N.rbc = @(u,v) [diff(v); v + 1];
[u, v] = N\0;
pass(6) = isa(u, 'chebfun') && isa(v, 'chebfun');

%% BVP two var case with info
dom  = [0, 2];
N = chebop(@(x,u,v) [diff(u,2) + 0.1*sin(u) + v; diff(v,2) + 0.1*cos(v) + u], dom);
N.lbc = @(u,v) [u - 1 ; diff(v)];
N.rbc = @(u,v) [diff(v); v + 1];
[u, v, info] = N\0;
pass(7) = isa(u, 'chebfun') && isa(v, 'chebfun') && isstruct(info);

%% IVP two var case
dom  = [0, 2];
N = chebop(@(x,u,v) [diff(u,2) + 0.1*sin(u) - v; diff(v,2) + 0.1*cos(v) - u], dom);
N.lbc = @(u,v) [u - 1 ; v + 1; diff(u); diff(v)];
u = N\0;
pass(8) = isa(u, 'chebmatrix');

%% IVP two var case, two outputs
dom  = [0, 2];
N = chebop(@(x,u,v) [diff(u,2) + 0.1*sin(u) - v; diff(v,2) + 0.1*cos(v) - u], dom);
N.lbc = @(u,v) [u - 1 ; v + 1; diff(u); diff(v)];
[u, v] = N\0;
pass(9) = isa(u, 'chebfun') && isa(v, 'chebfun');

%% IVP two var case with info
dom  = [0, 2];
N = chebop(@(x,u,v) [diff(u,2) + 0.1*sin(u) - v; diff(v,2) + 0.1*cos(v) - u], dom);
N.lbc = @(u,v) [u - 1 ; v + 1; diff(u); diff(v)];
[u, v, info] = N\0;
pass(10) = isa(u, 'chebfun') && isa(v, 'chebfun') && isstruct(info);

%% BVP four var case
dom  = [0, 2];
N = chebop(@(x,u,v,w,y) [ ...
    diff(u) + 0.1*sin(u) + w; ...
    diff(v) + 0.1*sin(w) + y;
    diff(w) + 0.1*sin(y) - u;
    diff(y) + 0.1*sin(u) - v], dom);
N.lbc = @(u,v,w,y) [u - 1 ; v + 1];
N.rbc = @(u,v,w,y) [w + .5; y - .5];
u = N\0;
pass(11) = isa(u, 'chebmatrix');

%% BVP four var case, multiple outputs
dom  = [0, 2];
N = chebop(@(x,u,v,w,y) [ ...
    diff(u) + 0.1*sin(u) + w; ...
    diff(v) + 0.1*sin(w) + y;
    diff(w) + 0.1*sin(y) - u;
    diff(y) + 0.1*sin(u) - v], dom);
N.lbc = @(u,v,w,y) [u - 1 ; v + 1];
N.rbc = @(u,v,w,y) [w + .5; y - .5];
[u,v,w,y] = N\0;
pass(12) = isa(u, 'chebfun') && isa(v, 'chebfun') && ...
    isa(w, 'chebfun') && isa(y, 'chebfun');

%% BVP four var case with info
dom  = [0, 2];
N = chebop(@(x,u,v,w,y) [ ...
    diff(u) + 0.1*sin(u) + w; ...
    diff(v) + 0.1*sin(w) + y;
    diff(w) + 0.1*sin(y) - u;
    diff(y) + 0.1*sin(u) - v], dom);
N.lbc = @(u,v,w,y) [u - 1 ; v + 1];
N.rbc = @(u,v,w,y) [w + .5; y - .5];
[u,v,w,y,info] = N\0;
pass(13) = isa(u, 'chebfun') && isa(v, 'chebfun') && ...
    isa(w, 'chebfun') && isa(y, 'chebfun') && isstruct(info);

%% IVP four var case
dom  = [0, 2];
N = chebop(@(x,u,v,w,y) [ ...
    diff(u) + 0.1*sin(u) + w; ...
    diff(v) + 0.1*sin(w) + y;
    diff(w) + 0.1*sin(y) - u;
    diff(y) + 0.1*sin(u) - v], dom);
N.lbc = @(u,v,w,y) [u - 1 ; v + 1; w-.5; y+.5];
u = N\0;
pass(14) = isa(u, 'chebmatrix');

%% IVP four var case, multiple outputs
dom  = [0, 2];
N = chebop(@(x,u,v,w,y) [ ...
    diff(u) + 0.1*sin(u) + w; ...
    diff(v) + 0.1*sin(w) + y;
    diff(w) + 0.1*sin(y) - u;
    diff(y) + 0.1*sin(u) - v], dom);
N.lbc = @(u,v,w,y) [u - 1 ; v + 1; w-.5; y+.5];
[u,v,w,y] = N\0;
pass(15) = isa(u, 'chebfun') && isa(v, 'chebfun') && ...
    isa(w, 'chebfun') && isa(y, 'chebfun');

%% IVP four var case with info
dom  = [0, 2];
N = chebop(@(x,u,v,w,y) [ ...
    diff(u) + 0.1*sin(u) + w; ...
    diff(v) + 0.1*sin(w) + y;
    diff(w) + 0.1*sin(y) - u;
    diff(y) + 0.1*sin(u) - v], dom);
N.lbc = @(u,v,w,y) [u - 1 ; v + 1; w-.5; y+.5];
[u,v,w,y] = N\0;
pass(16) = isa(u, 'chebfun') && isa(v, 'chebfun') && ...
    isa(w, 'chebfun') && isa(y, 'chebfun') && isstruct(info);

end
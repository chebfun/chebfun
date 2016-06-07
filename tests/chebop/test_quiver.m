function pass = test_quiver(~)
% TEST_QUIVER   Test the chebop/quiver command.
%
% chebop/quiver has to do with plotting, so we're mainly checking whether
% nothing crashes, not the correctness of the results.

hfig = figure('Visible','off');

%% Nonlinear pendulum, second order, default options
N = chebop(0, 10*pi);
N.op = @(t,y) diff(y, 2) + sin(y);
pass(1) = doesNotCrash(@() quiver(N,[-2 2 -1 1]));

%% Pass more options
pass(2) = doesNotCrash(@() quiver(N,[0 1 0 2],'xpts',25','ypts', 25, ...
    'scale',0.4,'normalize',true));

%% Lotka-Volterra, first order coupled system
N = chebop(@(t,u,v) [diff(u)-2.*u+u.*v; diff(v)+v-u.*v], [0 4]);
pass(3) = doesNotCrash(@() quiver(N, [0 2 0 4], 'normalize', true, ...
    'scale',.5,'linewidth',2));

%% Third order ODE should give an error
N = chebop(0, 10*pi);
N.op = @(t,y) diff(y, 3) + sin(y);
try 
    quiver(N,[-2 2 -1 1]);
catch ME
    pass(4) = strcmp(ME.identifier, 'CHEBFUN:CHEBOP:quiver:tooHighOrder');
end

%% Second order coupled system should also give an error
N = chebop(@(t,u,v) [diff(u,2)-2.*u+u.*v; diff(v)+v-u.*v], [0 4]);
try 
    quiver(N,[-2 2 -1 1]);
catch ME
    pass(5) = strcmp(ME.identifier, 'CHEBFUN:CHEBOP:quiver:tooHighOrder');
end

end


%% A function to test whether an expression crashes or not.
% Taken from tests/chebfun/test_plot.m
function pass = doesNotCrash(fn)

try
    fn();
    pass = true;
catch ME %#ok<NASGU>
    pass = false;
end

end

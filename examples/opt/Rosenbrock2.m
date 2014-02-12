%% The Rosenbrock function in 2D optimisation (revisited)
% Nick Hale, 4th March 2013

%%
% (Chebfun example opt/Rosenbrock2.m)
% [Tags: #Rosenbrock, #optimization, #Chebfun2]

LW = 'LineWidth';  lw = 1;
MS = 'MarkerSize'; ms = 20;
  
%% 1. The Rosenbrock function
% One of the very first Chebfun examples [2] demonstrated how standard 1D
% Chebfun often does a reasonable job of minimizing or maximizing a function
% defined on a 2D rectangle, such as the "Rosenbrock function" [1]:
f = @(x,y) (1-x).^2 + 100*(y-x.^2).^2;

%%
% With the introduction of Chebfun2, this task is made even more simple!
F = chebfun2(f, [-1.5 1.5 -1 3]);
[minf, minx] = min2(F)

%%
% We can plot the Rosenbrock function and the minimum we computed:
contour(F, 10:10:300, LW, lw), colorbar, shg
hold on, plot(minx(1), minx(2), '.k', MS, ms), hold off

%% 2. A function with several local minima
% The example [2] explained how Chebfun computed these results: by taking maxima
% along 1D slices, and then taking the maxiumum of these results. For functions
% with multiple local minima, like the one below, this meant the 'splitting'
% 'on' flag needed to be set.

f = @(x,y) exp(x-2*x.^2-y.^2).*sin(6*(x + y + x.*y.^2));

%% 
% Chebfun2 uses a different algorithm for locating maxima and minima, and so
% splitting is not required in this case:

F = chebfun2(f);
[minf, minx] = min2(F)

%%
% We can plot this function too, again with the computed minimum:
contour(F,30, LW, lw), colorbar, hold on
plot(minx(1), minx(2), '.k', MS, ms)

%%
% In fact, since Chebfun2 is really working with a 2D representation, we can
% investigate more deeply. 

%%
% For example we can compute all the stationary points of f and add these to our
% contour plot:
tp = roots(grad(F))
plot(tp(:,1), tp(:,2), 'kx', MS, 15, LW, 2)

%%
% We can make chebfuns of the zero level curves of F:
g = roots(F);
plot(g, 'k', LW, 2), shg

%%
% Or we can even make chebfuns of the .5 level curves of F:
g = roots(F - .5);
plot(g, 'r', LW, 2), shg

%%
% References:
%
% [1] H. H. Rosenbrock, "An automatic method for finding the greatest or least
% value of a function", Computer Journal, 3 (1960), 175-184.
%
% [2] http://www.chebfun.org/examples/opt/html/Rosenbrock.shtml

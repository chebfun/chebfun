%% The Dixon-Szego function in 2D optimisation
% Jari Fowkes and Nick Trefethen, November 2010

%%
% (Chebfun example opt/DixonSzego.m)
% [Tags: #optimization, #2D]
  
%%
% The Chebfun example opt/Rosenbrock.m
% shows how Chebfun can be used to minimize a function of two variables
% over a rectangle.  The present example is adapted from that one, and simply
% considers another function investigated by Dixon and Szego in 1975:
f = @(x,y) (4-2.1*x.^2+ x.^4/3).*x.^2 + x.*y + 4*(y.^2-1).*y.^2;

%%
% Over the rectangle [-2,2] x [-1.25,1.25], the function looks like this:
LW = 'linewidth'; FS = 'fontsize'; MS = 'markersize';
x = linspace(-2,2); y = linspace(-1.25,1.25);
[xx,yy] = meshgrid(x,y); ff = f(xx,yy);
figure, contour(x,y,ff,30,LW,1.2), colorbar
axis([-2 2 -1.25 1.25]), hold on

%%
% Here is Chebfun code taken from opt/Rosenbrock.m to find the minimum and
% plot the point where it is attained. (The minimum is actually achieved at
% two points, because of symmetry, but Chebfun does not detect this.)
tic
fminx0 = @(x0) min(chebfun(@(y) f(x0,y),[-1.25 1.25]));
fminx = chebfun(fminx0,[-2 2],'vectorize','splitting','on');
[minf,minx] = min(fminx)
[minf,miny] = min(chebfun(@(y) f(minx,y), [-1 3]))
toc
plot(minx,miny,'.k',MS,20)

%%
% Reference:
%
% [1] L. C. W. Dixon and G. P. Szego,
% The global optimization problem: an introduction, in
% L. C. W. Dixon and G. P. Szego (eds.),
% Towards Global Optimisation 2, North-Holland, Amsterdam 1978,
% pp. 1-15.

%% Restricted-denominator approximations
% Stefan Güttel, 18th April 2012

%%
% (Chebfun example approx/RestrictedDenominatorApproximations.m)
% [Tags: #rational, #exponential, #best]

%%
%
close all
clear all
LW = 'LineWidth'; lw = 2; MS = 'MarkerSize'; ms = 14;

%% 
% The approximation of the exponential function is of great practical
% interest in Scientific Computing. In a nutshell, this is based on the 
% fact that the solution of the ODE $y'(t) + Ay(t) = 0$ is given by 
% $y(t) = \exp(-tA)y(0),$ where $A$ is a square matrix and hence the 
% evaluation of $y(t)$ involves the computation of a matrix exponential.
% This is the basis for the design of so-called exponential integrators,
% which require efficient methods for computing the matrix exponential 
% (and other related functions). 
%
% If A is symmetric positive definite, it is reasonable to replace the 
% function $f(x) = \exp(-tx)$ by a rational approximation $r(x)$ on the 
% interval $[0,+\infty)$, and to evaluate $r(A)$ instead of $f(A)$. 

%%
% Let us first define the variable $x$, and the function $f$, setting the 
% time parameter $t$ equal to $1$ for the rest of this demo.

x = chebfun('x',[0,inf]);
f = exp(-x)

%%
% This chebfun f has a length of 41, and under the hood Chebfun has used
% a rational transformation to represent this function on an unbounded 
% domain: remember that every nonconstant polynomial would diverge to
% infinity on an unbounded domain, but f decays to 0 as x goes to infinity. 
% Mathematically, we have expanded f in rational Chebyshev functions [2].
% We can have a look at the associated map as follows

f.funs.map

%%
% This output tells us that the map is a rational transformation 
% $y = m(x) = (-15s + x)/(15s + x)$, mapping $[0,+\infty)$ to the interval
% $[-1,1)$. (The parameter a is the left endpoint of the unbounded domain, 
% i.e., a = 0 in our case.)  The parameter s is a scaling parameter that 
% can be changed. By default it is taken from mappref.parinf(1), in this 
% case s = 1. Therefore the rational transformation we are dealing with has 
% a pole at the point c = -15.  

%%
% Such rational approximations with prescribed poles often go under
% the name _restricted-denominator approximations_. One may ask what is the
% optimal pole c < 0 one must choose such that a restricted-denominator 
% rational function can achieve a smallest possible uniform approximation 
% error. To this end we consider an approximation of degree 4, with the 
% pole c varying in some negative interval. For each such parameter c we 
% compute the polynomial best uniform approximation to the exponential 
% function transplanted to [-1,1] using Chebfun's REMEZ command. 
% The following plot shows the uniform approximation error as a function 
% of c:

d = 4;
c = -logspace(-.5,1,100);
for j = 1:length(c), 
    f = chebfun(@(y) exp(c(j)*(y+1)./(1-y)),[-1,1]);
    p = remez(f,d);
    err(j) = norm(p - f,inf);
end
loglog(c,err,'b',LW,lw)
xlabel('parameter c')
ylabel('uniform approximation error')
axis([min(c),max(c),min(err)/3,max(err)])
legend('d = 4')
grid on

%%
% From this graph we can can read off that the best uniform approximation
% error of 3.14e-3 is obtained with the pole c = -5.7.
% In 2006, van den Eshof & Hochbruck conjectured that this graph
% has precisely d - 1 local minima, and even more, that these minima coincide
% with the local maxima of the graph for one degree higher [3]. We can also 
% verify this second part of their conjecture via Chebfun:

d = 5;
for j = 1:length(c), 
    f = chebfun(@(y) exp(c(j)*(y+1)./(1-y)),[-1,1]);
    p = remez(f,d);
    err(j) = norm(p - f,inf);
end
hold on
loglog(c,err,'r',LW,lw)
legend('d = 4','d = 5')

%%
% For d = 5 we read off an optimal pole c = -3.6, for which a uniform 
% approximation error of 1.09e-3 can be achieved. It has been proven by
% Andersson that the optimal concentrated pole c for approximating the
% exponential function converges to -inf when d is increased, and 
% behaves like $c \sim -d/\sqrt{2}$ [1]. 

%% 
% To gain some understanding what happens when two error graphs of 
% different order touch, let us investigate the error curve e = f - p 
% for d = 4 and a narrow parameter range, and then plot the points where
% the error curve attains its extremal values (the equioscillation points):

figure
d = 4;
c = -logspace(log10(2),log10(2.5),50);
for j = 1:length(c), 
    f = chebfun(@(y) exp(c(j)*(y+1)./(1-y)),[-1,1]);
    p = remez(f,d);
    e = f - p;
    [y,x] = minandmax(e,'local');
    ind = find(abs(y) > 0.99*max(abs(y)));
    x = x(ind);
    plot(x,c(j),'b.',MS,ms); hold on
end
xlabel('equioscillation points')
ylabel('parameter c')
axis([-1.2,1.2,min(c),max(c)])

%%
% Note that the d + 2 equioscillation points seem to have the tendency of 
% continuously moving leftwards as c goes to minus infinity. At the critical 
% parameter c = -2.23 the left-most point exits the interval [-1,1], and a 
% new point must enter from the right to keep the number of equioscillation 
% points at least d + 2. However, at this very moment of transition there 
% are d + 3 equioscillation points in the interval, which means that p also 
% is a polynomial best uniform approximation of degree d + 1, and therefore 
% the above graphs must meet for c = -2.23.

%%
% We would like to extend the conjecture by van den Eshof & Hochbruck [3] a 
% little further (and these observations were heavily stimulated by the
% ease with which Chebfun allows one to play around with functions):
% In fact, similar phenomena of touching error graphs can be observed
% for other polynomial approximations to the transplanted exponential 
% function. For example, here is the above plot with the uniform norm 
% replaced by the L2 norm:

figure
L = legpoly(0:5,'norm');
c = -logspace(-.5,1,100);
for j = 1:length(c), 
    f = chebfun(@(y) exp(c(j)*(y+1)./(1-y)),[-1,1]);
    p4 = L(:,1:5)*(L(:,1:5)'*f);
    err4(j) = norm(p4 - f,2);
    p5 = L(:,1:6)*(L(:,1:6)'*f);
    err5(j) = norm(p5 - f,2);
end
loglog(c,err4,'b',LW,lw); hold on
loglog(c,err5,'r',LW,lw)
xlabel('parameter c')
ylabel('L2 approximation error')
axis([min(c),max(c),min(err5)/1.5,max(err4)])
legend('d = 4','d = 5')
grid on

%%
% Again, each graph seems to have d - 1 local minima, and touches the graph
% of degree d at d - 1 points (but not necessarily at the minima or
% maxima). Similar graphs can be obtained for best weighted L2 approximants, 
% and even Chebyshev interpolants. A proof of these observations remains
% open, and one may wonder what special property of the (transplanted) 
% exponential function is required to cause these effects.

%%
% References:
%
% [1] J.-E. Andersson, Approximation of e^(-x) by rational functions with 
% concentrated negative poles, J. Approx. Theory (32), pp. 85--95, 1981.
% 
% [2] J. P. Boyd, Orthogonal rational functions on a semi-infinite
% interval, J. Comput. Phys. (70), pp. 63--88, 1987.
%
% [3] J. van den Eshof and M. Hochbruck, Preconditioning Lanczos 
% approximations to the matrix exponential, SIAM J. Sci. Comput. (27), 
% pp. 1438--1457, 2006.


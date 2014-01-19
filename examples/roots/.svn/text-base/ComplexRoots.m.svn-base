%% Roots of a complex function via Cauchy integrals
% Nick Trefethen, September 2011

%%
% (Chebfun example roots/ComplexRoots.m)
% [Tags: #roots, #complex, #Cauchy]

function ComplexRoots

%%
% Poles and zeros of complex functions can be located
% by the evaluation of contour integrals, as mentioned in
% Chapter 5 of the Chebfun Guide.  For example, 
% suppose we have a function like this one with a single
% root s1 in the unit disk:
ff = @(z) (z-0.5i).*exp(z);

%%
% We can find the root as the value of a contour integral
% around the unit circle:
%
%   s = (1/2i*pi) INT z (f'(z)/f(z)) dz

%%
% Since Chebfun works with real independent variables,
% we parametrize the unit circle by a real variable t on [-1,1]:
z = chebfun('exp(1i*pi*t)');

%%
% which gives us
%
%   s1 = (1/2i*pi) INT z ( (df/dt)(dt/dz) / f ) (dz/dt) dt
%
%      = (1/2i*pi) INT z ( (df/dt) / f ) dt

%%
% So here is the Chebfun evaluation:
f = ff(z);
s1 = sum(z.*diff(f)./f)/(2i*pi) 

%%
% There is nothing in this computation that depends on
% the use of the unit disk.  Other contours are equally tracatable
% in Chebfun, as illustrated in the Example complex/KeyholeContour and
% in Chapter 5 of the Chebfun Guide.

%%
% This method of finding a single root goes back at least
% to McCune in 1966 [3].  In practice we would often
% want to be able to find multiple roots, and a generalized algorithm
% for this case was published by Delves and Lyness in
% 1967 [1], with mathematical origins as far back as Jackson
% in 1917 [2].  The idea here is that if f has more
% than one root in the unit disk, then the value s above
% comes out as the sum of all these roots.  Similarly
%
%   s2 = (1/2i*pi) INT z^2 (f'(z)/f(z)) dz
%
% is the sum of the squares of the roots, the analogous
% formula for s3 with a factor z^3 gives the sum of the
% cubes, and so on.  And a count of the number of roots is
% given by
%
%   s0 = (1/2i*pi) INT (f'(z)/f(z)) dz
%
% (this is basically the argument principle).  So for example we can
% count the number of roots of cosh(pi*z) in the unit
% disk like this:

ff = @(z) cosh(pi*z);
f = ff(z);
s0 = sum(diff(f)./f)/(2i*pi)

%%
% Here are the sum of the roots and the sum of their squares:
s1 = sum(z.*diff(f)./f)/(2i*pi)
s2 = sum(z.^2.*diff(f)./f)/(2i*pi)

%%
% corresponding to roots at +- 0.5i.  We can find these numbers
% systematically by noting that the monic polynomial p(z)
% with these roots has coefficients c0=(s1^2-s2)/2, c1=-s1, c2=1.
% So here is a calculation of the two roots in the unit disk of cosh(pi*z):
p = [1 -s1 (s1^2-s2)/2];
roots(p)

%%
% Generalization to higher numbers of roots can be done
% via Newton's identities.  We don't pursue the general case
% here but instead
% write a code that finds three roots of an analytic function in
% the unit disk:

function r = roots3(ff)    % find 3 roots of ff in unit disk
z = chebfun('exp(1i*pi*t)');
f = ff(z);
s0 = sum(diff(f)./f)/(2i*pi);
s1 = sum(z.*diff(f)./f)/(2i*pi);
s2 = sum(z.^2.*diff(f)./f)/(2i*pi);
s3 = sum(z.^3.*diff(f)./f)/(2i*pi);
p = [1 -s1 (s1^2-s2)/2 -(s1^3-3*s1*s2+2*s3)/6];
r = roots(p);
end

%%
% Here is an example:
ff = @(z) cosh(exp(z)).*(z-.3).*(1+4*z.^2);
roots3(ff)

%%
% Here is another:
ff = @(z) (z.^3-1/8).*exp((-1-2i)*z);
roots3(ff)

end

%%
% References
%
% [1] L. M. Delves and J. N. Lyness, A numerical
% method for lcoating the zeros of an analytic function,
% Mathematics of Computation 21 (1967), 543-560.
%
% [2] D. Jackson, Roots and singular points of semi-analytic
% functions, Annals of Mathematics 19 (1917), 142-151.
%
% [3] J. E. McCune, Exact inversion of dispersion relations,
% Physics of Fluids 9 (1966), 2082-2084.

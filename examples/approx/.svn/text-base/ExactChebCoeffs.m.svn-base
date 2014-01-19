%% Exact Chebyshev expansion coefficients of a function
% Mark Richardson, 13th June 2012

%%
% (Chebfun example approx/ExactChebCoeffs.m)
% [Tags: #Chebyshevcoefficients, #residue]

%% 1. Introduction
% In this example, we shall compare the results of the Chebfun construction
% process to a known closed-form formula for the Chebyshev expansion
% coefficients of a function with a pole.
%
% This gives us an excellent excuse to visit some interesting approximation 
% theory from the 1960s!

%% 2. The residue method of Elliott
% For certain functions, explicit formulas for the coefficients in the
% Chebsyhev series expansion may be obtained using a contour-integral
% technique described by Elliott [1]. Here is how it works:
%
% The Chebsyhev coefficients of a Lipschitz-continuous function f can be 
% determined by the integral
% $$ a_n = \frac{2}{\pi} \int_{-1}^{1} \frac{f(x)T_n(x)}{\sqrt{1-x^2}} 
% {\rm d} x.$$ 

%%
% If f(x) is analytic within a particular contour C in the complex-plane, 
% then by Cauchy's integral formula, we can also write
% $$ f(x) = \frac{1}{2 \pi i} \int_{C} \frac{f(z)}{z-x} {\rm d} z.$$

%%
% If C is large enough to enclose the unit interval, then the second of 
% these two  formulas can be subsituted into the first and the orders of
% integration interchanged to give
%
% $$ a_n = \frac{1}{\pi^2 i} \int_C f(z) \int_{-1}^{1} 
% \frac{T_n(x) {\rm d} x}{(z-x)\sqrt{1-x^2}} {\rm d} z.$$

%%
% The integral with respect to x can be computed exactly so that we end up
% with
% $$a_n = \frac{1}{\pi i} \int_{C} \frac{f(z)}{\sqrt{z^2-1}(z \pm 
% \sqrt{z^2-1})^n} {\rm d} z.$$

%%
% Here, we note that $\rho = |z \pm \sqrt{z^2-1}|$ is the parameter of the 
% usual Bernstein Ellipse $E_\rho$ with the point $z = x + iy$ on its 
% boundary.
%
% So, if $f$ is a function with a pole at $z_0$ whose integral around 
% $E_\rho$ tends to zero as $\rho \to \infty$, then by the Residue theorem 
% we have
% $$ a_n = \frac{-2r_0}{\sqrt{z_0^2-1}(z_0 \pm \sqrt{z_0^2-1})^n},$$
% where $r_0$ is the residue of the pole at $z_0$.

%% 3. A function with  pole
% As an example, consider the function $$ f(x) = \frac{1}{5 + x} .$$

%%
% This function can be represented in Chebfun by an interpolant in 
% 17 points.
f  = @(x) 1./(5+x);
fc = chebfun(f);
k = 1:length(fc);
%%
% The function has a pole at -5 and corresponding residue 1. Substituing 
% these values into the above formula then gives the us following exact 
% expression for the Chebyshev expansion coefficients 
% $$ a_n =  \frac{1}{\sqrt{6}} \frac{(-1)^n}{(5+\sqrt{24})^n} $$

%%
% The theoretical coefficients match the ones computed by Chebfun, apart 
% from floating point representation and aliasing effects. The $a_0$ 
% coefficient is out by the usual factor of 2. 
exact_coeffs = flipud( (1/sqrt(6)*(-1).^(k-1)./(5+sqrt(24)).^(k-1))' );
cheb_coeffs = chebpoly(fc)';
display([exact_coeffs cheb_coeffs exact_coeffs-cheb_coeffs])

chebpolyplot(fc)
title('Chebyshev coefficients of 1/(5+x)')
xlabel('n'), ylabel('log(|a_n|)'), grid on
%%
% References:
%
% [1] D. Elliott, The evaluation and estimation of the coefficients in the
% Chebyshev series expansion of a function, Math. Comp. 18 (1964).

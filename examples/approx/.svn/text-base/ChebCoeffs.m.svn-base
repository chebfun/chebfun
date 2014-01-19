%% Chebyshev coefficients
% Nick Trefethen, September 2010

%%
% (Chebfun example approx/ChebCoeffs.m)
% [Tags: #Chebyshev, #coefficients, #CHEBPOLYPLOT]

%%
% Every function defined on [-1,1], so long as it is at least Lipschitz
% continuous, has an absolutely and uniformly convergent Chebyshev series:
%
%     f(x) = a_0 + a_1 T_1(x) + a_2 T_2(x) + ....
%
% The same holds on an interval [a,b] with appropriately scaled and shifted
% Chebyshev polynomials.

%%
% For many functions you can compute these coefficients with the command
% CHEBPOLY.  For example, here we compute the Chebyshev coefficients of a
% cubic polynomial:
x = chebfun('x');
format long
disp('Cheb coeffs of 99x^2 + x^3:')
p = 99*x.^2 + x.^3;
a = chebpoly(p)'

%%
% Notice that following the usual Matlab convention, the coefficients
% appear in order from highest degree to lowest.  Thus it is often more
% useful to write
disp('Cheb coeffs of 99x^2 + x^3:')
a = chebpoly(p)'; a = a(end:-1:1)

%%
% Similarly, here are the Chebyshev coefficients down to level 1e-15 of
% exp(x):
disp('Cheb coeffs of exp(x):')
a = chebpoly(exp(x))'; a = a(end:-1:1)

%%
% You can plot the absolute values of these numbers on a log scale with
% CHEBPOLYPLOT:
FS = 'fontsize'; MS = 'markersize'; LW = 'linewidth';
chebpolyplot(exp(x),'.-',LW,1,MS,16), grid on
xlabel('degree n',FS,12)
ylabel('|a_n|',FS,12)
title('Chebyshev coefficients of exp(x)',FS,16)

%%
% Here's a similar plot for a function that needs thousands of terms to be
% represented to 15 digits.  (Can you explain why it looks like a wide
% stripe?)
chebpolyplot(exp(x)./(1+10000*x.^2)), grid on
xlabel('degree n',FS,12)
ylabel('|a_n|',FS,12)
title('Chebyshev coefficients of exp(x)/(1+10000x^2)',FS,16)

%%
% These methods will work for any function f that's represented by a global
% polynomial, i.e., a chebfun consisting of one fun. (Normally this means
% that the Chebyshev series needs fewer than 65536 terms, Chebfun's default
% value of MAXDEGREE. To change this default, type HELP CHEBFUNPREF.) What
% about Chebyshev coefficients for functions that are not smooth enough for
% such a representation?  Here one can use the TRUNC option in the Chebfun
% constructor. For example, suppose we are interested in the function
f = sign(x);
figure, plot(f,'k',LW,2), ylim([-1.5 1.5])
title('sign(x)',FS,16)

%%
% If we try to compute all the Chebyshev coefficients, we'll get an error.
% On the other hand we can compute the first ten of them like this:
p = chebfun(f,'trunc',10);
a = chebpoly(p)'; a = a(end:-1:1)

%%
% Here's the degree 9 polynomial obtained by adding up these first terms of
% the Chebyshev expansion:
hold on
plot(p,'m',LW,2)
title('sign(x) and truncated Chebyshev series',FS,16)

%%
% This is not the same as the degree 9 polynomial interpolant through 10
% Chebyshev points:
pinterp = chebfun(f,10);
plot(pinterp,'--','color',[0 .8 0],LW,2)
title('Same, also with Chebyshev interpolant',FS,16)

%%
% Reference
%
% [1] L. N. Trefethen, Approximation Theory and Approximation Practice,
% draft book available at http://www.maths.ox.ac.uk/chebfun/ATAP/.

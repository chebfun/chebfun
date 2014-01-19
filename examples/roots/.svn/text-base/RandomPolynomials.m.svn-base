%% Roots of random polynomials
% Nick Trefethen, 6th November 2011

%%
% (Chebfun example roots/RandomPolynomials.m)
% [Tags: #roots, #random]

%%
% If p(z) = a_0 + a_1 z + ... + a_n z^n is a polynomial with random
% coefficients, its roots tend to lie near the unit circle [1,2].
% To be specific, suppose a_0,...,a_(n-1) are independent real
% numbers from the standard normal distribution, with a_n=1.
% Here are typical pictures for n = 50 and 200, based on Matlab's
% ROOTS command.  Incidentally, such computations are numerically
% stable [4];
% the well-known difficulties of polynomial rootfinding pertain mainly
% to problems with roots away from the unit circle and coefficients
% far from random.
nn = [50 200];
MS = 'markersize'; FS = 'fontsize'; ms = 12; fs = 12;
for j = 1:2
    subplot(1,2,j)
    n = nn(j);
    a = [1; randn(n,1)];
    r = roots(a);
    plot(r,'.k',MS,ms)
    axis(1.5*[-1 1 -1 1]), axis square
    title(['monomial, n=' int2str(n)],FS,fs)
end

%%
% The reason these roots lie near the unit circle is that
% the monomials 1, z, z^2,... are orthogonal polynomials defined on
% that set.  For orthogonal polynomials defined on a different region
% of the complex plane, we get roots tending to cluster on the
% boundary of that region [3].  In particular, the best-known families
% of orthogonal polynomials are defined over [-1,1], and random
% polynomials expressed in these bases tend to have roots near
% [-1,1].  The Chebfun ROOTS command provides an easy way to compute
% roots of such polynomials stably.
% For example, here is what happens if we repeat the 
% experiment above but for random polynomials in the basis of
% Chebyshev polynomials, i.e., p = a_0 + a_1 T_1 + ... + a_n T_n.
for j = 1:2
    subplot(1,2,j)
    n = nn(j);
    a = [1; randn(n,1)];
    p = chebfun(a,'coeffs');
    r = roots(p,'all','norecurse');
    plot(r,'.k',MS,ms)
    axis(1.5*[-1 1 -1 1]), axis square
    title(['Chebyshev, n=' int2str(n)],FS,fs)
end

%%
% Legendre polynomials are orthogonal over the same interval [-1,1],
% so the results aren't much different:
for j = 1:2
    subplot(1,2,j)
    n = nn(j);
    a = [1; randn(n,1)];
    A = legpoly(n:-1:0);
    p = A*a;
    r = roots(p,'all','norecurse');
    plot(r,'.k',MS,ms)
    axis(1.5*[-1 1 -1 1]), axis square
    title(['Legendre, n=' int2str(n)],FS,fs)
end

%%
% References:
%
% [1] J. Hammersley, The zeros of a random polynomial,
% Proceedings of the Third Berkeley Symposium on Mathematical
% Statistics and Probability, 1954-1955, U. California Press,
% 1956, pp. 89-111.
%
% [2] L. A. Shepp and R. J. Vanderbei, The complex zeros of
% random polynomials, Transactions of the American Mathematical
% Society 347 (1995), 4365-4384.
%
% [3] B. Shiffman and S. Zelditch, Equilibrium distribution of
% zeros of random polynomials, International Mathematical Research
% Notices, 2003, pp. 25-49.
%
% [4] K.-C. Toh and L. N. Trefethen, Pseudozeros of polynomials and
% pseudospectra of companion matrices, Numerische Mathematik, 68
% (1994), 403-425.

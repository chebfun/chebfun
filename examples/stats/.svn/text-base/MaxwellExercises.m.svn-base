%% Maxwell Distribution -- Exercises from a Textbook
% Jie Gao, 23rd September 2013

%%
% (Chebfun example stats/MaxwellExercises.m)
% [Tags: #Maxwell distribution, #probability, #mean, #variance]

%% 1. Introduction
% Probability and statistics textbooks contain many exercise problems concerning
% various probability distributions.  In this Example we use Chebfun to solve a
% problem involving the Maxwell distribution from the textbook [1]. We use a
% different form of the probability density function from the one given in the
% textbook to make it easier to compare the expectation with the exact value.
% Maxwell distributions are often observed and used in physical applications,
% such as describing particle speeds in gases for a certain temperature. Other
% similar Examples look at problems from the same book involving the normal,
% exponential, uniform, gamma, beta, and Rayleigh distributions.

%%
% Like most textbooks, [1] emphasizes problems that can be solved on paper and
% don't need numerical tools such as Chebfun. As soon as one varies the problem
% a little, however, numerical solutions often become necessary. Thus first we
% solve the problem as written, and then we examine an application inspired by
% [2].

%% 2. Problem 20, page 127

%%
% _The Maxwell distribution is defined by
% $$f(x;b)=\frac{\sqrt{2}}{b^3\sqrt{\pi}}x^2\exp(\frac{-x^2}{2b^2}).$$ Calculate
% the mean and variance._

%%
% Let us take b to be 2.3.
b = 2.3;
%%
% The domain is the positive real line because the probability of X being less
% than zero is zero under a Maxwell distribution. The probability density
% function (PDF) of X can be defined like this:
ff = @(x) sqrt(2)/(b^3*sqrt(pi)).*x.^2.*exp(-x.^2/(2*b^2));
f = chebfun(ff, [0, inf]);

%%
% We can compute the expectation of X like this:
format long
x = chebfun(@(x) x, [0, inf]);
meanm = sum(chebfun(@(x) x.*ff(x), [0, inf]))
%%
% We compare this value to the exact mean value:
meanm_exact = 2*sqrt(2/pi)*b
%%
% The expectation of X.^2 is:
meanmsqr = sum(chebfun( @(x) x.^2.*ff(x), [0, inf]))
%%
% So the variance is:
varm = meanmsqr - meanm^2
%%
% The standard deviation is thus:
stdm = sqrt(varm)
%%
% Let us plot the PDF, the mean, and the one-standard-deviation region of this
% Maxwell distribution:
hold off, plot(f, 'k', 'linewidth', 1.6, 'interval', [0, 10]), 
hold on, x = (meanm - stdm):.1:(meanm + stdm);
h = area(x,f(x));
set(h, 'FaceColor', [0.1 0.6 0.7])
plot( [meanm, meanm], [0, f(meanm)], '--y' )
grid on

%% 3. Problem 20, page 127 -- application
%%
% As we have said earlier, the Maxwell distribution can describe the
% distribution of the speed of a molecule. Since this is a continuous
% distribution, we cannot obtain the probability of a molecule having a speed of
% exactly 3.0 m/s, but we can determine the probability that a particle's speed
% is between 2.9 m/s and 3.1 m/s by integrating the PDF using these two values
% as the limits.
%%
% Let us take b = 1.3:
b = 3.7;
%%
% As before, we have the PDF of the speed of a particle to be:
ff = @(x) sqrt(2)/(b^3*sqrt(pi)).*x.^2.*exp(-x.^2/(2*b^2));
f = chebfun(ff, [0, inf]);
%%
% The cumulative density function (CDF) of X is the indefinite integral of
% f:
fint = cumsum(f);
%%
% Now we calculate the desired probability:
P = fint(3.1) - fint(2.9)
%%
% We can evaluate the average speed of this particle:
v_avg = sum(chebfun(@(x) x.*ff(x), [0, inf]))
%%
% Let us plot the PDF, the mean, and the one-standard-deviation region of 
% this Maxwell distribution:
hold off, plot(f, 'k', 'linewidth', 1.6, 'interval', [0, 18]), 
hold on, x = 2.9:.1:3.1;
h = area(x,f(x));
set(h, 'FaceColor', [0.1 0.6 0.3])
plot( [v_avg, v_avg], [0, f(v_avg)], '--r' )
grid on
%%
% The average speed-squared is:
vsqr_avg = sum(chebfun(@(x) x.^2.*ff(x), [0, inf]))
%%
% Now we can compute the average kinetic energy of a particle of mass 0.15 kg by
% the formula K = 1/2*m*v^2:
m = 0.15;
K_avg = 1/2*0.15*vsqr_avg
%%
% Reference:
% 
% [1] A. M. Mood, F. A. Graybill, and D. Boes, Introduction to the Theory of
% Statistics, McGraw-Hill, 1974.
%
% [2] John Baliga, Maxwell-Boltzmann Distribution. Available at:
% http://academics.triton.edu/faculty/jbaliga/pdfs/MaxwellBoltzmannDistribution.pdf
% (Accessed: 15 Sep 2013)
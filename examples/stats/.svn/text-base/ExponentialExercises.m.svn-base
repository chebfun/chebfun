%% Exponential Distribution -- Exercises from a Textbook
% Jie Gao and Nick Trefethen, 6th May 2013

%%
% (Chebfun example stats/ExponentialExercises.m)
% [Tags: #exponentialdistribution, #probability]

%% 1. Introduction
% Probability and statistics textbooks contain many exercise problems concerning
% various probability distributions.  In this Example we use Chebfun to solve
% two problems involving the exponential distribution from the textbook [1].
% The exponential distribution describes the lifetimes of things of many
% different kinds. It models the length of the time interval between successive
% events occurring in a Poisson process. Other similar Examples look at problems
% from the same book involving the normal, beta, gamma, Rayleigh, and Maxwell
% distributions.

%%
% Like most textbooks, [1] emphasizes problems that can be solved on paper and
% don't need numerical tools such as Chebfun. As soon as one varies the problem
% a little, however, numerical solutions often become necessary. Thus first we
% solve the problems as written, and then we solve a variant for each problem.
% For the standard exponential distribution, we also show an application adapted
% from [2].

%% 2. Problem 1(g), page 124

%%
% _If a random variable X has a negative exponential distribution with mean 2,
% find P[X<1|X<2]._

%%
% The probability density function (PDF) of the exponential distribution can be
% defined like this:
ff = @(x,lambda) lambda.*exp(-lambda.*x);

%%
% Since the mean of the distribution is 2, we know the value of lambda is 1/2.
lambda = 1/2;
%%
% We choose the domain to be the interval [0, inf] since the probability of the
% random variable X being less than zero is zero under an exponential
% distribution. So we can construct the chebfun like this:
f = chebfun(@(x) ff(x,lambda), [0, inf]);

%%
% The cumulative density function (CDF) is the indefinite integral of f:
fint = cumsum(f);

%%
% By Bayes' Theorem, P[X<1|X<2] = P[X<1]/P[X<2]. Hence we can find the
% probability P[X<1|X<2] like this:
format long
p = fint(1)/fint(2)

%%
% Here is a graph of the normalized f on the interval [0, inf] and the region
% X<1.
hold off, g = f/fint(2);
h = area(g{0, 1});
set(h,'FaceColor',[0.2 0.25 0.6]), axis auto
LW = 'linewidth';
hold on, plot(g,'k',LW,1.6,'interval',[0 2]), grid on

%% 3. Problem 1(g), page 124 -- numerical variant

%%
% Before we move on to a numerical variant, let us first compute P[X>8|X>5] and
% compare it to P[X>3] to give an example of showing an important property of
% exponential distribution -- memorylessness. That is, P[X>s+t|X>s] = P[X>t] for
% s, t >0.

%%
% By Bayes' Theorem, P[X>8|X>5] = P[X>8]/P[X>5]. We know that P[X>a] = 1-P[X<=a]
% for all a >= 0. So this conditional probability can be computed like this:
cp = (1-fint(8))/(1-fint(5))

%%
% We now compare the value of q to P[X>3] and we shall see that they are
% essentially equal:
q = 1-fint(3)

%%
% For the numerical variant, let us now change the linear term in the previous
% PDF to a log(2) power and compute the probability P[X<1|X<3]:
f = chebfun(@(x) exp(-lambda.*x.^(log(2))), [0, inf], 'splitting', 'on');

%%
% We normalize the function to make it a proper density function:
f = f/sum(f);

%%
% The desired probability is computed similarly as before:
fint = cumsum(f);
p = fint(1)/fint(3)

%%
% Now let us make a plot showing the probability:

hold off, f = f/fint(3);
h = area(f{0, 1});
set(h,'FaceColor',[0.35 0.9 0.4]), axis auto
LW = 'linewidth';
hold on, plot(f,'b',LW,1.6,'interval',[0 3]), grid on

%% 4. Problem 1(g), page 124 -- application

%%
% _Suppose the lifespan in hundreds of hours, T, of a light bulb of a home lamp
% is exponentially distributed with lambda = 0.2, i.e. the PDF of T is f(t) =
% 0.2*exp(-0.2*t) for t>0. lambda = 0.2 is called the failure rate of the light
% bulb. The reliability at time t is defined by Rel(t) = P[T>t]._

%%
% We first define t as a chebfun:
t = chebfun('t', [0, inf]);

%%
% So the PDF of T is defined like this:
f = 0.2*exp(-0.2*t);
%%
% The CDF is the indefinite integral of f:
fint = cumsum(f);
%%
% We can define the reliability like this:
Rel = 1-fint(t);
%%
% Let us compute the probability that the light bulb will last more than 700
% hours:
p1 = 1-fint(7)
%%
% Also, the probability that the light bulb will last more than 900 hours is:
p2 = 1-fint(9)
%%
% We now plot Rel to examine the reliability of the light bulb as t increases:
hold off, plot(Rel, 'k', 'linewidth', 1.6), axis auto, grid on

%%
% Suppose we want to know when the reliability is 10%. That is, we need to
% choose d such that P[X>d] = 0.1:
d = roots(1-fint-0.1)
%%
% We can see that there is an approximately 90% likelihood that the light bulb
% will last less than 1150 hours. Equivalently, 10% of the time the light bulb
% will last more than 1150 hours.

%%
% Let us plot this probability:
hold off, h = area(f{d, inf});
set(h,'FaceColor',[0.5 0.1 0.1]), axis auto
LW = 'linewidth';
hold on, plot(f,'k',LW,1.6), grid on

%% 5. Problem 1(l), page 124

%% 
% _Suppose X has a negative exponential distribution with parameter lambda. If
% P[X<=1] = P[X>1], what is var[X]?_

%% 
% This is equivalent to P[X<=1] = 1-P[X<=1], i.e. P[X<=1] = 1/2. In other
% words, median = 1.

%%
% The probability density function (PDF) of the exponential distribution can be
% defined like this:
ff = @(x,lambda) lambda.*exp(-lambda.*x);

%%
% Similarly, we define the PDF on the positive real line and construct the
% chebfun like this:
f = @(lambda) chebfun( @(x) ff(x, lambda), [0, inf]);

%%
% We can calculate lambda which satisfies P[X<=1] = P[X>1] like this:
lambda = roots( chebfun( @(lambda) sum(f(lambda),[0,1])-.5, [0, 10],'vectorize'))
%%
% We compare this value to the exact value of lambda, log(2):
lambda_exact = log(2)
%%
% We can find the expectation of X like this:
meanx = sum(chebfun( @(x) x.*lambda.*exp(-lambda.*x), [0,inf]))
%%
% Correspondingly, the expectation of X^2 is:
meanxsqr = sum(chebfun( @(x) x.^2*lambda.*exp(-lambda.*x), [0,inf]))
%%
% Thus, the variance of the distribution is:
var = meanxsqr - meanx^2
%%
% We can compare this value to the exact value of the variance, 1/(log(2))^2:
var_exact = 1/(log(2))^2
%%
% Hence, the standard deviation of the distribution is:
stdexp = sqrt(var)
%%
% Let fn be the new function with lambda known:
fn = chebfun(@(x) lambda*exp(-lambda.*x),[0, inf]);
%%
% We now plot fn and the standard deviation of the distribution:
hold off, plot( [stdexp, stdexp], [0, 1], '--r', 'linewidth', 2),
hold on, plot(fn, 'k', 'linewidth', 1.6, 'interval', [0, 8]), grid on

%% 6. Problem 1(l), page 124 -- numerical variant

%%
% Let us first do a similar computation, except replacing the linear term in the
% previous exponential distrubtion by an absolute value with a 13/5 power times
% log(x+1/2). Because of the absolute value and the zero of log(x+1/2) at x =
% 1/2, this distribution has a singularity at x = 1/2. We take lambda to be 1.5
% and construct the chebfun like this:
lambda = 1.5;
f = chebfun( @(x) exp(-abs((lambda.*x.^(13/5).*log(x+1/2)))),...
[0, inf], 'splitting', 'on');

%%
% We normalize the distribution as before:
f = f/sum(f);

%%
% We can compute the variance of X in the following steps:

%%
% The expectation of X is:
meanx = sum(chebfun(@(x) x.*f(x), [0, inf], 'splitting', 'on'))
%%
% The expectation of X.^2 is:
meanxsqr = sum(chebfun( @(x) x.^2.*f(x), [0,inf], 'splitting', 'on'))
%%
% Now we have the variance:
var = meanxsqr - meanx^2
%%
% And the standard deviation is thus:
stdexp = sqrt(var)

%%
% We can also compute the mode:
[fmax,mode] = max(f)
%%
% Now we plot the PDF of the distribution, its standard deviation (in cyan
% dashed line), mean (in magenta -. line), and mode (marked by two blue stars).
% Note the singularity at x = 1/2:
hold off, plot( [stdexp, stdexp], [0, 1], '--c', 'linewidth', 2),
hold on, plot ([meanx, meanx], [0, 1], '-.m', 'linewidth', 2),
plot ([mode, mode], [0, 1], 'p'),
plot(f, 'b', 'linewidth', 1.6, 'interval', [0, 2]), grid on

%%
% References:
% 
% [1] A. M. Mood, F. A. Graybill, and D. Boes, Introduction to the Theory of
% Statistics, McGraw-Hill, 1974.
%
% [2] Prof. Jane M. Horgan, Chapter 17 Applications of the Exponential
% Distribution. Available at: www.computing.dcu.ie/~jhorgan/chapter17slides.pdf
% (Accessed: 6 May 2013)

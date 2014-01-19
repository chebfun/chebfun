%% Beta Distribution -- Exercises from a Textbook
% Jie Gao, 2nd July 2013

%%
% (Chebfun example stats/BetaExercise.m)
% [Tags: #Beta distribution, #probability, #mode, #median,
% #Bayesian inference, #Hypothesis Testing, #prior odds, #posterior odds]

%% 1. Introduction
% Probability and statistics textbooks contain many exercise problems concerning
% various probability distributions.  In this Example we use Chebfun to solve a
% problem involving the beta distribution from the textbook [1].  Beta
% distribution is a family of probability densities taking on values in the
% interval (0, 1). It can model experiments of a variety of shapes. Other
% similar Examples look at problems from the same book involving the normal,
% exponential, uniform, gamma, Rayleigh, and Maxwell distributions.

%%
% Like most textbooks, [1] emphasizes problems that can be solved on paper and
% don't need numerical tools such as Chebfun. As soon as one varies the problem
% a little, however, numerical solutions often become necessary. Thus first we
% solve the problem as written, and then we solve a variant. Lastly we
% demostrate how the beta distribution can be useful in Bayesian inference by an
% example from [3].


%% 2. Problem 2(a), page 124

%%
% _Find the mode of the beta distribution._

%%
% We know that the probability density function (PDF) of a beta distribution can
% be defined like this: f = B.*x^(a-1).*(1-x)^(b-1), where constant = 1/h(a,b),
% a>0, b>0; and B is Beta function: the definite integral of
% x^(a-1).*(1-x)^(b-1) with respect to X from 0 to 1. a and b are two parameters
% controlling the shape of the distribution.

%%
% As the probability of X being less than zero or bigger than one is zero under
% a beta distribution, let us first make a chebfun of X on the interval [0, 1]:
x = chebfun('x', [0, 1]);

%%
% Let us take a to be 1.5 and b to be 3. We can construct the numerator of the
% PDF like this:
a = 1.5;
b = 3;
f = x.^(a-1).*(1-x).^(b-1);
%%
% Let the denominator of the PDF, the constant B, be:
B = sum(f);
%%
% The probability density function (PDF) of the beta distribution is thus:
f = f/B;
%%
% Since the mode is simply the value of X when f(X) reaches its maximum, we use
% the Chebfun command MAX to find the mode:
format long
[fmax,mode] = max(f)
%%
% We can compare this value to the exact value of the mode:
mode_exact = (a-1)/(a+b-2)
%%
% Let us plot the PDF and the mode of this beta distribution:
LW = 'LineWidth';
hold off, plot(f, 'k', LW, 1.6), axis auto
hold on, plot( [mode,mode], [0, 2], '--r', LW, 1.6 ), grid on
%%
% We can see that the density function is skewed, since the two parameters that
% determine the shape of the density function (i.e. a and b) are not equal to
% each other.
%% 3. Problem 2(a), page 124 -- numerical variant

%%
% In this problem, we want to calculate the mode and median of the following
% distribution. We choose the domain to be [0, 2]. Let us replace the beta
% distribution by changing X to X/2, 1-X to 2-X/2, and the powers of these two
% terms from a-1 and b-1 to log(a+1) and exp(b-1):
g = chebfun(@(x) (1/2*x).^(log(a)).*(2-1/2*x).^(exp(b-1).*sqrt(x)),...
[0, 2], 'splitting', 'on');
%%
% As before, we need to normalize the function by hand:
g = g/sum(g);
%%
% The mode is calculated like this:
[gmax,gmode] = max(g);
display(gmode)
%%
% Since we want to find the median of this distribution, we need to know the
% value of m in P[X<m] = 1/2. We use the cumulative distribution function (CDF)
% to solve for m:
gint = cumsum(g);
m = roots(chebfun(@(m) gint(m)-0.5, [0, 2]))
%%
% Let us plot the PDF, the mode (the green dashed line), and the median (the
% black dashed line) of this distribution:
hold off, plot(g, 'b', LW, 1.6), axis auto
hold on, plot( [gmode,gmode], [0, 1], '--g', LW, 1.6 ), 
plot( [m, m], [0, 1], '--k', LW, 1.6), grid on
%%
% We can see that the mode and median are very close to each other. In fact, the
% graph of the probability density function of this distribution is rather
% similiar to a normal distribution.

%% 4. Problem 2(a), page 124 -- application in Bayesian inference

%%
% In Bayesian inference, we treat unknown parameters (e.g. theta) as random
% variables. A prior density summarizes information about this unknown parameter
% theta before observing the data x. We can write the prior density as h(theta).
% Beta distributions can be prior probability distributions for geometric,
% binomial, and Bernoulli distributions. Correspondingly, a posterior density -
% a conditional density of theta given x - contains updated information about
% theta and can be written as h(theta|x). We can write the probability model for
% data x as f(x|theta) (i.e., the likelihood function) since we want to
% emphasize its conditionality on theta [1].
%%
% By Bayes' Theorem, an important relation among the prior, posterior, and
% likelihood function is: posterior is proportional to the product of likelihood
% and prior, i.e. h(theta|x) is proportional to f(x|theta) * h(theta).
%%
% In many situations, we want to compare two hypotheses H_0 and H_1, exactly one
% of which is true. The following example shows how Bayesian inference acts in
% hypothesis testing.

%%
%
%%
%
%%
% _Example from [3]_
%%
% We have two products: Product P_0 is old and standard; Product P_1 is newer
% and more expensive.
%%
% We make the following assumptions: (1) the probability theta that a customer
% prefers P_1 has prior h(theta) which is Beta(a,b) (2) the number of customers
% X (out of n) that prefer P_1 is X ~ Binomial(n,theta).
%%
% Let's say theta >= 0.6 means that P_1 is a substantial improvement over P_0.
% So take H_0: theta >= 0.6 and H_1: theta < 0.6.
%%
% We consider 3 possible priors: Beta(0,5,0.5) -- Jeffreys' prior Beta(1,1) --
% uniform prior Beta(2,2) -- "skeptical" prior (i.e. favors values of theta near
% 1/2). Let us plot these three prior densities in Chebfun:
%% 
% We can define a chebfun in terms of theta:
theta = chebfun('theta', [0, 1]);
%%
% We first define the Jeffreys' prior density, f_j:
a_j = 0.5;
b_j = 0.5;
f_j = theta.^(a_j-1).*(1-theta).^(b_j-1);
B_j = sum(f_j);
f_j = f_j/B_j;
%%
% We then define the uniform prior density, f_u:
a_u = 1;
b_u = 1;
f_u = theta.^(a_u-1).*(1-theta).^(b_u-1);
B_u = sum(f_u);
f_u = f_u/B_u;
%% 
% We next define the skeptical prior density, f_s:
a_s = 2;
b_s = 2;
f_s = theta.^(a_s-1).*(1-theta).^(b_s-1);
B_s = sum(f_s);
f_s = f_s/B_s;
%%
% Now we can plot them in one graph:
hold off, plot(f_j, 'k', f_u, '--r', f_s, '--g', LW, 1.6), 
ylim([0 3.5]), grid on
hleg1 = legend('Beta(0.5,0.5) Jeffreys prior', ...
    'Beta(1,1) uniform prior',...
    'Beta(2,2) skeptical prior');
xlabel('theta')
ylabel('prior density')
%%
% Prior odds = P(H_0)/P(H_1), where P(H_0) = definite integral of probability
% density function of the corresponding beta distribution for theta from 0.6 to
% 1, and P(H_1) = definite integral of probability density function of the
% corresponding beta distribution for theta from 0 to 0.6. [The odds of any
% event A are P(A)/(1-P(A)).]
%%
% So we can define the cumulative density function (CDF) for each density function
% as usual:
f_jint = cumsum(f_j);
f_uint = cumsum(f_u);
f_sint = cumsum(f_s);
%%
% Now we compute P(H_1) for each of the three priors:
P_H_1_j = f_jint(0.6);
P_H_1_u = f_uint(0.6);
P_H_1_s = f_sint(0.6);
%%
% Then we compute P(H_0) for each prior:
P_H_0_j = 1 - P_H_1_j;
P_H_0_u = 1 - P_H_1_u;
P_H_0_s = 1 - P_H_1_s;
%%
% So the prior odds of each prior are:
prior_odds_j = P_H_0_j/P_H_1_j
prior_odds_u = P_H_0_u/P_H_1_u
prior_odds_s = P_H_0_s/P_H_1_s
%%
% Suppose we have x = 13 "successes" from n = 16 customers. Then the posterior
% h(theta|x) is Beta(x+a,n-x+b) (we call this 'new') with x = 13 and n = 16.
% (Note: a can be a_j, a_u, or a_s; b can be b_j, b_u, or b_s.) To avoid
% confusion, we use y instead of x here:
y = 13;
n = 16;
%%
% Posterior odds = P(H_0|x)/P(H_1|x), where P(H_0|x) = definite integral of
% probability density function of the corresponding new beta distribution for
% theta from 0.6 to 1, and P(H_1) = definite integral of probability density
% function of the corresponding new beta distribution for theta from 0 to 0.6.
%%
% We first define the new beta distribution for the Jeffreys' prior:
f_j_n = theta.^(y+a_j-1).*(1-theta).^(n-y+b_j-1);
B_j_n = sum(f_j_n);
f_j_n = f_j_n/B_j_n;
%%
% We then define the new beta distribution for the uniform prior:
f_u_n = theta.^(y+a_u-1).*(1-theta).^(n-y+b_u-1);
B_u_n = sum(f_u_n);
f_u_n = f_u_n/B_u_n;
%%
% We next define the new beta distribution for the skeptical prior:
f_s_n = theta.^(y+a_s-1).*(1-theta).^(n-y+b_s-1);
B_s_n = sum(f_s_n);
f_s_n = f_s_n/B_s_n;
%%
% We define the cumulative density function (CDF) for each density function as
% before:
f_j_nint = cumsum(f_j_n);
f_u_nint = cumsum(f_u_n);
f_s_nint = cumsum(f_s_n);
%%
% Then we compute P(H_1|x) for each of the three priors (Note: To avoid
% confusion in Matlab, we write 'given' instead of the symbol '|' in our
% statements):
P_H_1_j_n = f_j_nint(0.6);
P_H_1_u_n = f_u_nint(0.6);
P_H_1_s_n = f_s_nint(0.6);
%%
% Also we compute P(H_0|x) for each prior:
P_H_0_j_n = 1 - P_H_1_j_n;
P_H_0_u_n = 1 - P_H_1_u_n;
P_H_0_s_n = 1 - P_H_1_s_n;
%%
% So the posterior odds of each prior is:
pos_odds_j = P_H_0_j_n/P_H_1_j_n
pos_odds_u = P_H_0_u_n/P_H_1_u_n
pos_odds_s = P_H_0_s_n/P_H_1_s_n
%%
% We plot the corresponding new beta distributions for the Jeffreys', uniform,
% and skeptical priors, respectively:
hold off, plot(f_j_n, 'k', f_u_n, '--r', f_s_n, '--g', LW, 1.6), 
ylim([0 4.5]), grid on
hleg2 = legend('Beta(13.5,3.5)', 'Beta(14,4)', 'Beta(15,5)');
set(hleg2, 'Position', [.3,.5,.1,.2]);
xlabel('theta')
ylabel('posterior density')
%%
% The Bayes factor is the ratio of the posterior odds to the prior odds.
bay_f_j = pos_odds_j/prior_odds_j
bay_f_u = pos_odds_u/prior_odds_u
bay_f_s = pos_odds_s/prior_odds_s
%%
% Let us make a table for three types of priors, their prior odds, posterior
% odds, and Bayes factors.
f = figure('Position',[100 100 400 150]);
% I can explicitly enter the values of the above three variables and use
% uitable:
dat = {'Beta(0.5,0.5)', 0.7728, 26.6073, 34.4317;
    'Beta(1,1)', 0.6667, 20.5411, 30.8116;
    'Beta(2,2)', 0.5432, 13.3650, 24.6037;};
cnames = {'Prior', 'Prior odds', 'Posterior odds', 'Bayes factor'};
columnformat = {'char', 'numeric', 'numeric', 'numeric'};
t = uitable('Position', [20 20 360 100], 'Data', dat,...
    'ColumnName', cnames,'ColumnFormat', columnformat, 'RowName', []);
%%
% If we want to to call the variable values directly from the computation we can
% use a different command to make the table - fprintf:
betaArgs = [0.5, 0.5; 1.0, 1.0; 2.0, 2.0];
col1 = [prior_odds_j; prior_odds_u; prior_odds_s];
col2 = [pos_odds_j; pos_odds_u; pos_odds_s];
col3 = [bay_f_j; bay_f_u; bay_f_s];
fprintf('==================================================================\n');
fprintf( '%s           %s      %s      %s\n', 'Prior', 'Prior odds', ...
    'Posterior odds', 'Bayes factor');
fprintf('==================================================================\n');
nRows = size(col1,1);
nCols = 3;
for i = 1:nRows
    fprintf( 'Beta(%.1f,%.1f)    %.6f        %.6f           %.6f\n', ...
    betaArgs(i,1), betaArgs(i,2), col1(i), col2(i), col3(i)); 
end
fprintf('==================================================================\n');
    

%%
% By [2], we can conclude that we have strong evidence for H_0 since all three
% Bayes factors are between 20 and 150.

%% 
% Acknowledgements:
% 
% I am thankful to Prof. Nick Trefethen whose suggestions for improving the code
% have been very useful. I also thank Mohsin Javed for allowing me to use his
% code of fprintf.

%% 
% References:
% 
% [1] A. M. Mood, F. A. Graybill, and D. Boes, Introduction to the Theory of
% Statistics, McGraw-Hill, 1974.
%
% [2] Neil Laws, Part A Statistics: Bayesian Inference, Hilary Term 2013.
% Available at: http://www.stats.ox.ac.uk/~laws/partA-stats/A-StatsBayes.pdf
% (Accessed: 1 Jul 2013)
%
% [3] Neil Laws, Example from Carlin and Louis (2008), Bayesian Inference.
% Available at:
% http://www.stats.ox.ac.uk/~laws/partA-stats/2013A-BayesSlides.pdf (Accessed: 1
% Jul 2013)
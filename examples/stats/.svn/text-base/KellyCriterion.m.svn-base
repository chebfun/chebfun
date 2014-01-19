%% Optimal bet-sizing and the Kelly Criterion
% Mark Richardson, 24th October 2012

%%
% (Chebfun example statistics/KellyCriterion.m)
%
% [Tags: #KellyCriterion, #Optimisation, #Investment]
%%
LW = 'LineWidth'; MS = 'MarkerSize'; format short

%% 1. Theoretical setup 
% 
% Suppose you have a fixed amount of money to invest, say £100, and 
% you are made aware of a particular investment opportunity which may 
% be entered into an unlimited number of times, each time on identical 
% terms. The rules of the wager are as follows. Your total balance grows 
% or shrinks over time depending on the accumulative outcomes of your bets. 
% You may bet any amount of your balance on each wager. With probability 
% $p$, your stake is trebled and returned to you; with probability $1-p$, 
% your stake is lost. Two questions follow: 
% 
% 1) Do you want to play this game?
% 
% 2) If so, how much do you want to wager?
% 
% Regarding the first of these questions: Notwithstanding any philosophical 
% or other objections, the answer is, clearly, that it depends on $p$. 
% In fact, we can see quite easily that since the payoff odds 
% of this  wager sit at $2:1$ (i.e., you can win two units for every one 
% unit wagered), then if $p \leq 1/3$ then we would certainly not wish to 
% play since the expected value of the game is nonpositive. So let us 
% assume that $p > 1/3$ in order that the answer to the first question is 
% "yes". In fact, for definiteness, let's say $p = 1/2$.
%
% An answer to the second question was first given in 1956 by a Bell Labs 
% scientist named John Kelly [1].  (Kelly, incidentally was working at the 
% time on problems in the fledgling discipline of information theory 
% together with eminent scientists such as Claude Shannon of the famous 
% Nyquist-Shannon sampling theorem.) Kelly's solution to the bet-sizing 
% problem became known as the "Kelly Criterion". It was used to great 
% effect later on by, amongst others, the MIT Blackjack card-counting team. 
% It is now a mainstay of any self-respecting first course on basic 
% investment principles. 
% 
% To set up this problem, we consider a quantity $G$ called the 
% "logarithmic growth rate of capital". Suppose that our starting amount 
% is $C_0$, and that $C_n$ is the amount we have accumulated after the 
% $n^{th}$ play of the game. Then, we define
%
% $$ G := \lim_{n \to \infty} \frac{1}{n} \log \left( \frac{C_n}{C_0} \right) .$$
%
% Each time, we shall seek to wager a fixed proportion $0 \leq f \leq 1$ 
% of our capital. Suppose that $w$ is the number of wins and $l$ is the 
% number of losses resulting from $n$ wagers. Note that $n = w + l$. 
% Suppose also that our game is generalised such that the payoff odds of 
% the game are $a:1$ for some $a > 0$, so that we win $a$ units with 
% probability $p$, and lose $1$ unit with probability $1-p$. Then we have
%
% $$ C_n = (1+af)^w (1-f)^l C_0 .$$
%
% Thus,
%
% $$ G = \lim_{n \to \infty}\left( \frac{w}{n}\log(1+af) + 
%       \frac{l}{n}\log(1-f) \right), $$
%
% so that, with probability one, we obtain
%
% $$ G(f) = p\log(1+af) + (1-p)\log(1-f).$$
% 
% This is the logarithm of the expected increase in our capital per bet. 
% Our task now is to attempt to maximise this quantity over $f \in [0,1]$. 
% Working through some basic calculus (i.e., differentiating and computing 
% the root of the resulting linear function), we find that $G$ is maximised 
% at  
%
% $$ f = \frac{ap - (1-p)}{a} .$$ 
%
% This is the Kelly Criterion. For our problem above, we had $a = 2$,
% $p = 1/2$. Using the Kelly Criterion formula in this case gives us 
% $f = 1/4$. In other words, the optimal strategy is to wager one quarter 
% of our bankroll each time we play.
% 
% Let's look at the expected capital growth function $\exp(G(f))$ for this 
% example. Note that we use a truncated interval due to the logarithmic 
% singularity at $x=1$. 
p = 0.5; a = 2;
G = @(f) p*log(1+a*f) + (1-p)*log(1-f);
GG = chebfun(G,[0 0.999]);
plot(exp(GG),LW,2), grid on
[Gmax,f] = max(GG); 
hold on, plot(f,exp(Gmax),'.r',MS,20), hold off
title('Expected rate of capital growth per bet')
xlabel('f'), ylabel('exp(G(f))'), ylim([0 1.2])
%%
% The computed optimal fraction agrees with the Kelly formula:
disp(f)
%%
% Using this strategy over the long term, we can therefore expect, on 
% average, our capital to be multiplied by the following quantity per bet:
disp(exp(Gmax))
%% 
% With these game parameters, betting anything less than a quarter of our
% wealth is too conservative, and betting any more is too aggressive. 
% We also note additionally that we expect to make money over the long 
% term so long as $G(f)$ is positive. The crucial interval for which this 
% is the case can be trivially determined using roots:
roots(GG)

%% 2. Numerical approach 
%
% Numerical methods were not really need above; they just provided a useful
% way for us to check the theory. However, should we wish to attack more 
% complicated problems, numerical methods, and in particular Chebfun, will 
% be of value. Indeed, the binomial problem considered above is 
% really far too simplistic to be of much use in practical applications. 
%
% The following example involving staggered payoffs is still artifical, 
% but is hopefully closer to a situation one might encounter in 
% practice.
% 
% We have a bet with six possible outcomes. Associated with each 
% outcomes are six probabilies $p_j$, which sum to 1, and some associated 
% payoffs $a_j$. Suppose that the probabilities and payoffs are as follows:
p1 = 0.4; p2 = 0.21; p3 = 0.26; p4 = 0.1; p5 = 0.02; p6 = 0.01;
a1 = -1 ; a2 = 0;    a3 = 1.2;  a4 = 1.3; a5 = 1.4;  a6 = 10;
%%
% We can then follow the same setup as before by defining:
G = @(f) p1*log(1+a1*f) + p2*log(1+a2*f) + p3*log(1+a3*f) + ...
            p4*log(1+a4*f) + p5*log(1+a5*f) + p6*log(1+a6*f);
%%
% Now the global optimisation step is taken care of numerically. This
% gives us the following:
GG = chebfun(G,[0 0.999]);
plot(exp(GG),LW,2), grid on
[Gmax,f] = max(GG); 
hold on, plot(f,exp(Gmax),'.r',MS,20), hold off
title('Expected rate of capital growth per bet')
xlabel('f'), ylabel('exp(G(f))'), ylim([0 1.2])
%%
% And the correct proportion to bet in this case is therefore:
disp(f)
%% 
% Finally, we should expect to lose our wealth if we bet any more than the 
% following fraction:
r = roots(GG); disp(r(2))

%% 3. Other comments
% 
% We examined this problem in the somewhat idealised setting for which
% accurate and a priori knowledge of the probabilities and corresponding 
% payoffs was assumed. However, such information is rarely available in 
% practice; typically proababilities are estimated from past data. For 
% this reason, it is common to take a more conservative approach and bet 
% less than the Kelly fraction -- "half-Kelly" is typical. This approach 
% also has the additional benefit of reducing the volatility (read 
% annualised standard-deviation of logarithmic returns) which can result 
% from betting the full Kelly fraction.
%
%%
% References:
%
% [1] J. L. Kelly, Jr, A New Interpretation of Information Rate, Bell 
%       Systems Technical Journal, 35, (1956), 917–-926
% 
% [2] Unknown author: http://www.elem.com/~btilly/kelly-criterion/
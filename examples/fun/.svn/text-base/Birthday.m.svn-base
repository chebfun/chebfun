%% Birthday cards and analytic functions
% Nick Trefethen, September 2010

%%
% (Chebfun example fun/Birthday.m)
% [Tags: #complex, #gift]

%%
% Chebfun's SCRIBBLE command was introduced for
% entertainment, but it turns out to be surprisingly useful
% also illustrating complex variables.  Suppose for example
% it is Chebyshev's birthday and you want to send him a card:
s = scribble('Happy Birthday Pafnuty!');
LW = 'linewidth'; lw = 1.8;
plot(s,'-',LW,lw)
xlim([-1.1 1.1]), axis equal, grid on

%%
% This chebfun s is a a piecewise linear
% complex funtion of a real variable, as we can see
% by writing it without the semicolon:
s

%%
% Since s is a chebfun, we can apply functions to it.
% For example, here is exp(s):

plot(exp(s),'b',LW,lw), axis equal, grid on

%%
% Here is exp(3i*s):
plot(exp(3i*s),'m',LW,lw), axis equal, grid on

%%
% Playing around with different functions is a good way
% to learn about complex variables, and a good way to
% make greeting cards.  Here are a couple more with axes
% turned off for greater beauty.
plot(exp((1+1i)*s),'g',LW,lw), axis equal, axis off
snapnow
plot(sinh(3*s),'r',LW,lw), axis equal, axis off


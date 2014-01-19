%% Delta and Heaviside Hyperfunctions
% Mohsin Javed, 17th April 2012

%%
% (Chebfun example complex/Hyperfuns.m)
% [Tags: #DIRAC, #complex, #deltafunction, #heaviside, #hyperfunction, #delta]

%% Introduction
% Hyperfunction theory describes a generalized function f(x) on the real
% line by viewing it as the difference between two generating functions
% F+(z) and F-(z), holomorphic in the upper and lower complex plane,
% respectively. The value of f(x) at a "regular point" x on the real line
% is calculated by taking the difference of the upper and lower part of the
% generating functions and then applying the limit:
%%
% limit ep-->0 [ F+(x+i*ep) - F-(x-i*ep) ]
%%
% i.e. we approach the real-axis from directly above and directly below the
% point x. If the limit does not exist, x is called a "singular point" of f
% and it does not make sense to talk about the value of f at x [1].

%%
% (A given hyperfunction does not have a unique choice of generators; you
% can add any analytic function to F+ and F- without changing their
% difference along the real axis.  So properly speaking, a hyperfunction is
% an equivalence class.)

%% The Dirac-Delta Function
% An elegant choice of generating function for the Dirac-delta function is to
% take both F- and F+ equal to the same function,
F = @(z) -1./(2i*pi*z);

%%
% Thus in this case F+ = F-.

%% The Heaviside Function
% The Heaviside function also has an elegant choice of generator, again
% with G- = G+ equal to the same function, the integral of the previous
% one:
G = @(z) -1/(2*pi*1i)*log(-z);

%%
% Using the above definitions, we define anonymous functions to describe
% the delta function and the Heaviside function on the interval [-1,1].
% Note that we take the real part only to remove imaginary rounding errors.
x = chebfun('x');
hyperDelta = @(ep) real(F(x+1i*ep)-F(x-1i*ep));
hyperHeaviside = @(ep) real(G(x+1i*ep)-G(x-1i*ep));


%%
% We now imitate the limiting process by taking small values of the
% parameter ep and plot the corresponding functions below:

%%
%
for ep = .1:-.01:.001;
    subplot( 1,2,1)
    plot(hyperDelta(ep)), hold on
    subplot( 1,2,2)
    plot(hyperHeaviside(ep)), hold on
end
subplot(1,2,1), title( 'Delta function' );
subplot(1,2,2), title( 'Heaviside function' );

%%
% Reference:
%
% [1] Urs Graf, Introduction to Hyperfunctions and Their Integral
% Transforms: An Applied and Computational Approach, Birkhaeuser, 2010.

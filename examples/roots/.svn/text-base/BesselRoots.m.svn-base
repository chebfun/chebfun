%% Roots of a Bessel function
% Nick Trefethen, September 2010

%%
% (Chebfun example roots/BesselRoots.m)
% [Tags: #roots, #Bessel]

%%
% Here is the Bessel function J0 on the interval [0,100].
J0 = chebfun(@(x) besselj(0,x),[0 100]);
figure, plot(J0,'linewidth',1.6), grid on
title('Bessel function J_0','fontsize',16)

%%
% We can find its roots like this:
r = roots(J0);
hold on, plot(r,J0(r),'.r','markersize',20)

%%
% The number of roots can be found with the LENGTH command:
number_of_roots = length(r)

%%
% Suppose you wanted to know the numbers of roots in
% various intervals [a,b].
% You could define an anonymous function:
rootsab = @(a,b) length(roots(chebfun(@(x) besselj(0,x),[a b])));

%%
% For example:
tic
disp('Number of roots between 1000000 and 1001000:')
n = rootsab(1000000,1001000)
toc



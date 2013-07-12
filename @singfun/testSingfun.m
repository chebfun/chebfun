%%
clear classes

a = 1.5; b = 1.0;
op = @(x) sin(1-x)./((1-x).^(b).*(1+x).^(a));
exponents = [-a, -b];
%exponents = [];
isSingEnd = [1,1];
singType = {'branch', 'branch'};
pref = [];
f = singfun( op, exponents, isSingEnd, singType, pref )
%%
imagf = imag(f)

%% 
realf = real(f)
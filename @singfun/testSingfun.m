%%
clear classes

a = .5; b = 1.5;
op = @(x) sin(1-x)./((1-x).^(b).*(1+x).^(a));
exponents = [-a, -b];
%exponents = [];
isSingEnd = [1,1];
singType = {'branch', 'branch'};
pref = [];
f = singfun( op, exponents, isSingEnd, singType, pref )
plot(f)
%%
imagf = imag(f)

%% 
realf = real(f)

%%
f = singfun( @(x) sin(2*pi*(x+1))./(x+1), [], [1 0], {'pole', 'pole'}, []  )
plotData(f)


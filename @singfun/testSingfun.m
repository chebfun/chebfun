%%
clear classes

a = .5; b = 1.5;
op = @(x) sin(1-x)./((1-x).^(b).*(1+x).^(a));
exponents = [-a, -b];
%exponents = [];
singType = {'sing', 'sing'};
pref = [];
f = singfun( op, exponents, singType, pref )

%%
imagf = imag(f)

%% 
realf = real(f)

%%
f = singfun( @(x) sin(2*pi*(x+1))./((1-x).^3.*(x+1).^2), [], {'pole', 'pole'}, []  )

%%
f = singfun( @(x) sin(x)./((1-x).^3.5.*(1+x).^.5), [], {'sing', 'sing'}, [] );
g = singfun( @(x) cos(x)./(1-x), [], {'sing', 'sing'}, [] );
s = f + g
xx = -.9:.01:.9;
error = feval(s, xx) - (feval(f, xx)+ feval( g, xx));
norm(error, inf )
%plot(xx, error)

%%
f = singfun(@(x) 1./(1-x), [], {'sing', 'sing'})
fp = diff(f);
xx = -.99:.01:.99;
error = feval(fp, xx) - 1./(1-xx).^2;
norm(error, inf )
%plot(xx, error)

%%
a = pi;
f = singfun(@(x) sin(a*x)./(1-x).^2, [], [], [] );
fp = diff(f);
xx = -.99:.01:.99;
fpExact = @(x) a*cos(a*x)./(1-x).^2 + 2*sin(a*x)./(1-x).^3;
error = feval(fp, xx) - fpExact(xx);
norm(error, inf )
%plot(xx, error)
%%
% works but if we change f to -f, it doesn't
f = chebfun( @(x) sin(12*pi*(1-x))./(1-x).^2, 'blowup', 'on' )
plot(f) % the plot is wrong?
[a, b] = min(f)
%%
% f changed to -f
f = chebfun( @(x) -sin(1-x)./(1-x).^2, 'blowup', 'on' )
a = min(f)

%%
f = singfun( @(x) sin(12*pi*(1-x))./(1-x).^2 )
g = 1./f
%% 
f = singfun( @(x) (1-x).*(1+x) )
1./f
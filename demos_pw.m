clear classes
d = [-1 0.5 1];
I = linop.eye(d);
D = linop.diff(d);
Z = linop.zeros(d);
x = chebfun('x',d);

%%
% One var with breakpoint
L = chebmatrix({D^2+I});
L = lbc(L,0);
L = rbc(L,1);
f = chebfun(1,d);
u = L\f
plot(u{1})
shouldbeSmall = norm( L{1,1}*u{1} - f )

%%
% Two vars with breakpoint
L = [ D^2-I, D^2+I; D^2+I, I-D ];
L = lbc(L,[I Z],0);
L = lbc(L,[Z I],-1);
L = rbc(L,[I Z],0);
L = rbc(L,[Z I],1);
f = chebmatrix({abs(x-0.5);sin(x)});
u = L\f
clf
plot(u{1}), hold on, plot(u{2},'k')
r = L*u - f;
norm(r{1})

%%
% Test automatic inclusion of breakpoints from RHS data.
d = [-1 1];
I = linop.eye(d);
D = linop.diff(d);
Z = linop.zeros(d);
x = chebfun('x',d);
f = chebmatrix({abs(x)});
L = chebmatrix({D^2-I});
L = lbc(L,0);
L = rbc(L,D,1);
u = L\f

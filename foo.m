function foo

problem = 0;

% coupled system of odes
if ( problem == 0)
dom = [-1,1];
x = chebfun('x', dom);

N = chebop(dom);
N.op = @(x,u1,u2,u3) [ diff(u1,2) + exp(x).*u3;...
                       diff(u2,2) + x.^2.*u1;...
                       diff(u3,2) + sqrt(2+x).*u2 ];
N.rbc = @(u1,u2,u3)  [diff(u1);u2;diff(u3)]; 
N.lbc = @(u1,u2,u3)  [u1;diff(u2);u3]; 
%N.bc = 'periodic';

Nstar = adjoint(N);
nargin(N.op)

f = sin(x);
g = cos(30*x);
h = g./(1.1+f);
rhs = [f;g;h];
%rhs = [ sin(pi*x); cos(2*pi*x); exp(sin(pi*x)+cos(pi*x)) ];
subplot(3,1,1), plot(rhs)

u = N\rhs;
subplot(3,1,2), plot(u)

v = Nstar\rhs;
subplot(3,1,3), plot(v)

ip1 = v'*rhs
ip2 = rhs'*u
commutator = abs(ip1 - ip2)/max(norm(v),norm(u))

% scalar ivp control problem
elseif ( problem == 1)
dom = [0,10];
a = 1;
b = 1;
alpha = 1;
beta = 0;
x = chebfun('x', dom);

g = @(u,w) .5*(a*u.^2 + b*w.^2);
gu = @(u,w) a*u;
gw = @(u,w) b*w;

f = chebop(@(x,u,w) diff(u,2) + sin(u) - w, dom);
f.lbc = @(u,w) [u-alpha; diff(u)-beta]; 

u0 = alpha + beta*(x-dom(1));
u0 = 0*x;
w0 = 0*x;

%[u,w,gval] = fmincon(g,f,u0,w0,gu,gw);
[u,w,gval] = fmincon(g,f,u0,w0);

% scalar bvp control problem
elseif ( problem == 2 )
    
dom = [-1,1];
a = 1;
b = 1;
alpha = 1;
beta = 1.3;
x = chebfun('x', dom);

g = @(u,w) .5*(a*u.^2 + b*w.^2);
gu = @(u,w) a*u;
gw = @(u,w) b*w;

f = chebop(@(x,u,w) .01*diff(u,2) + cos(u) - w, dom);
f.lbc = @(u,w) u-alpha; 
f.rbc = @(u,w) u-beta; 

u0 = 0*x;
w0 = 0*x;

[u,w,gval] = fmincon(g,f,u0,w0);

end


% u = u0; 
% w = w0;
% for ii = 3:5
% 
% n = 2^ii+1;
% [x,w] = chebpts(n,dom);
% y = [u(x);w(x)];
% 
% G = @(y) sum(g(chebfun(y(1:n),dom),chebfun(y(n+1:2*n),dom)));
% 
% nonlcon = @(y) cons(y,f,n);
% 
% options = optimoptions('fmincon','Display','iter','Algorithm','sqp');
% y = fmincon(G,y,[],[],[],[],[],[],nonlcon,options);
% 
% u = chebfun(y(1:n),dom);
% w = chebfun(y(n+1:2*n),dom);
% subplot(2,1,1), plot([u,w])
% subplot(2,1,2), plotcoeffs([u,w]), pause(5)
% 
% end
% 
% function [c,ceq] = cons(y,f,n)
% dom = f.domain;
% c = [];
% ceq = [ feval(f([chebfun(y(1:n),dom);chebfun(y(n+1:2*n),dom)]),x);...
%         feval(f.lbc(chebfun(y(1:n),dom),chebfun(y(n+1:2*n),dom)),dom(1))];
% end



%n = 64;
%cheboppref.setDefaults('minDimension',n,'maxDimension',n);
%cheboppref.setDefaults('bvpTol',1e-6);



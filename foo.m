function foo

problem = 0;

% coupled system of odes
if ( problem == 0)
dom = [-1,1];
x = chebfun('x', dom);

N = chebop(dom);
N.op = @(x,u1,u2) [ diff(u1,2) + x.*diff(u2) + u2; u1 + (2+cos(x)).*diff(u2,2) ];
%N.op = @(x,u1) [ diff(u1,2); cos(x).*u1 + diff(u1,2) ];
%N.op = @(x,u1,u2) [ diff(u1,2) + u2; u1 + diff(u2,2) ];
%N.op = @(x,u1,u2) [ diff(u1,2)+diff(u2); diff(u2,2) ];
N.lbc = @(u1,u2)  [diff(u1);u2]; 
N.rbc = @(u1,u2)  [u1;diff(u2)+u2]; 

L = linearize(N);
[Lstar,op] = adjoint(L,'bvp');

f = -sin(10*x);
g = cos(30*x);
rhs = [f;g];

pref = cheboppref();
pref.discretization = @chebcolloc2;
u = linsolve(L,rhs,pref);
u1 = u{1}; u2 = u{2};
subplot(2,1,1), plot([u1,u2])

v = linsolve(Lstar,rhs,pref);
v1 = v{1}; v2 = v{2};
subplot(2,1,2), plot([v1,v2])

v'*(L*u)
(Lstar*v)'*u
commutator = abs(v'*(L*u) - (Lstar*v)'*u)/max(norm(v),norm(u))

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



clc

A = chebop(@(x,u) diff(u,2) + sin(2*pi*x)*u + sum(u).*u);
A.bc = 'periodic'; 
u = A\1;
norm(A*u-1)

f = chebfun(@(x) sin(5*pi*x));
K = @(x,y) cos(-pi*(x-y));
A = chebop(@(x,u) 1/25*diff(u,2) + fred(K, u));
A.bc = 'periodic'; 
u = A\f;
norm(A*u-f)

%%

f = chebfun(@(x) sin(5*pi*x));
K = @(x,y) cos(-pi*(x-y));
A = chebop(@(x,u) 1/25*diff(u,2) + volt(K, u));
A.bc = 'periodic'; 
u = A\f;
norm(A*u-f)

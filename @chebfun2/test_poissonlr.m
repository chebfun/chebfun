%example/test for low rank poisson code. 
%%
% test 1: zero dirichlet condition
% rank 1 RHS: use fadi:
   f = chebfun2( @(x,y) exp(x +y) );

% for now set discretization: 
N = 128; 
%%
% solve with poisson: 
ut = chebfun2.poisson(f, 0, N,N);
   
%%
% solve with poisson_lr: 

u = chebfun2.poisson_lr(f, 0, N, N); % 0 = use fadi, not fiadi.


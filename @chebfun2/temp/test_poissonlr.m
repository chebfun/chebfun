%example/test for low rank poisson code. 
%%
% test 1: zero dirichlet condition (basic accuracy)

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

%%
% test 2: rank 1 RHS timings: 
K = 10; 
N = round(logspace(2.3, 3.7, K)); 
ranks = [1, 10, 50, 100, 200]; 
t_fullrank = zeros(K,1); 
t_lowrank = zeros(K,length(ranks)); 
%%
for j = 1:K
    t1 = tic; 
    ut = chebfun2.poisson(f, 0, N(j),N(j));
    tt = toc(t1); 
    t_fullrank(j) = tt; 
    j
end
%%
for j = 1:length(ranks)
    f = chebfun2(f, rank(j));
    for k = 1:K
    t2 = tic; 
    u = chebfun2.poisson_lr(f, 0, N(k), N(k)); % 0 = use fadi, not fiadi.
    tt = toc(t2); 
    t_lowrank(k,j) = tt; 
    [k, j]
    end
end

%%
% picture
loglog(N.', t_fullrank, 'k--', 'Linewidth', 2.5)
hold on
loglog(N.', t_lowrank.', 'Linewidth', 2.5)
legend('adi', 'rank = 1', 'rank = 10', 'rank = 50', 'rank = 100', 'rank =200', ...
    'Location', 'Northwest')
title('Discretization size vs. timings')
set(gca, 'fontsize', 18)
set(gcf, 'color', 'white')
hold off
shg
    
    








%% Absolute value approximations by rationals 2
% Yuji Nakatsukasa,  31st July 2012

%%
% (Chebfun example approx/AbsoluteValueScaled.m)
% [Tags: #rational, #Newton, #ABS]

%%
% This is a follow-up of Example approx/AbsoluteValue [1]. The goal is to
% find a rational approximation to the absolute value function $|x|$. That
% example used  Newton's method applied to $x^2=r^2$ with the initial guess
% $r=1$, given by the iteration $r := (r^2+x^2)/2r$. After $k$ steps we
% have a rational function of type $(2^k,2^k)$, which approaches $|x|$ as
% $k\rightarrow \infty$.

%%
% Let's rerun the code in that example and plot the error:

LW = 'linewidth'; lw = 1.6; FS = 'fontsize'; fs = 12;
x = chebfun('x',[-1 0 1]);
r = chebfun('1',[-1 0 1]);
kmax = 5; % # of iterations
for k = 0:kmax
    r = (r.^2+x.^2)./(2*r);
end
%%
semilogy(abs(r-abs(x))+eps,LW,lw), grid on
axis([-1 1 1e-18 10]),xlabel('x',FS,fs)
title('Error',FS,fs)

%%
% The main issue here is that the error is large near the origin, given
% that  the optimal type $(2^k,2^k)$ rational approximants to $|x|$
% achieves root-exponential accuracy $O(\mbox{exp}(-C\sqrt{2^k}))$ in the
% infinity norm [5,6].

%%
% Here we try another approach, which is to combine the formula
% $|x|=x/\mbox{sign}(x)$ with the scaled Newton iteration for approximating
% the sign function $\mbox{sign}(x)$. Newton's iteration for
% $\mbox{sign}(x)$ is defined by $r := (r+1/r)/2$ and the scaled Newton
% iteration is its scaled variant $r := (tr+1/(tr))/2$, where $t>0$ is
% determined so as to optimize the convergence. It requires a parameter
% $0<b<1$ such that the sign function is approximated on the interval
% $[b,1]$. For details on the scaled Newton iteration for the sign function
% see for example [2],[3, Ch. 8]. Once $r$ approximates $\mbox{sign}(x)$
% well, we get an approximation to $|x|$ via $r:=x/r$. As above, after $k$
% steps we have a type $(2^k,2^k)$ rational function that approximates
% $|x|$. Let's see how it works with an example:

rs = chebfun('x',[-1 0 1]);b=1e-3;
t=1/sqrt(b);
for k = 0:kmax
    if k>0
        t=sqrt(2/(t+1/t)); % scaling t
    end
    rs= ((t*rs)+1./(t*rs))/2;
end % rs now approximates the sign function
rs=x./rs; % get approximant to abs(x) via abs(x)=x/sign(x)
hold on, semilogy(abs(rs-abs(x)),'r',LW,lw)
axis([-1 1 1e-18 10]), grid on
xlabel('x',FS,fs)
title('Error',FS,fs)
legend('Newton','scaled Newton','location','best')

%% 
% Now the error is uniformly small across the interval $[-1,1]$. In fact,
% it can be shown that for a given $k$, the scaled Newton iteration yields
% the type $(2^k,2^k-1)$ best rational approximation to $\mbox{sign}(x)$ on
% the interval $[b,1]$ due to Zolotarev. Since the best type $(n,n)$
% approximation to $\mbox{sign}(x)$ yields accuracy
% $O(\mbox{exp}(-C\sqrt{n}))$ [5, Ch.4], we can show that also for $|x|$,
% the above process (with an appropriately chosen $b$) yields the optimal
% accuracy $O(\mbox{exp}(-C\sqrt{2^k}))$.
% 

%%
% The asymmetry, also observed in the example approx/AbsoluteValue [1],
% seems more pronounced in the above red plot. This is due to rounding
% errors: to observe this, let's see the plots for varying $k$.

clf;
semilogy(abs(r-abs(x)),LW,lw),hold on, % plot unscaled Newton
b=1e-3;kmax=5;
r = x; % initialize
colork=['k','m','g','r'];
t=1/sqrt(b);
for k = 0:kmax
    if k>0
        t=sqrt(2/(t+1/t)); 
    end
    r= ((t*r)+1./(t*r))/2;    
    if 1<k 
    semilogy(abs(x./r-abs(x)),'Color',colork(k-1),LW,lw),hold on        
    end
end 
legend('Newton k=5','s-Newton k=2','s-Newton k=3',...
    's-Newton k=4','s-Newton k=5','location','best')
axis([-1 1 1e-18 10]), grid on

%%
% Clearly for $k\leq 3$ the error is symmetric about the imaginary axis,
% exhibiting a near-equioscillating property. It is still curious that in
% the red plot, the effect of rounding error is present at a much larger
% value than the machine precision $10^{-16}$.


%%
% References:
%
% [1] http://www.chebfun.org/examples/approx/html/AbsoluteValue.shtml
%
% [2] R. Byers and H. Xu. A new scaling for Newton's iteration for the
% polar decomposition and its backward stability. SIAM J. Matrix Anal.
% Appl., 30(2):822-843, 2008.
%
% [3] N. J. Higham. Functions of Matrices: Theory and Computation. SIAM,
% Philadelphia, PA, USA, 2008.
%
% [4] D. J. Newman, Rational approximation of abs(x), Michigan Mathematical
% Journal 11 (1964), 11-14.
%
% [5] P. P. Petrushev and V. A. Popov, Rational Approximation of Real
% Functions, Cambridge University Press, 2011.
%
% [6] L. N. Trefethen, Approximation Theory and Approximation Practice,
% SIAM, to appear in late 2012.




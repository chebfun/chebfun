function [p,xi,sumnorms,sumshistory] = L1min_watson(f,n,xiini,reps)
% attempts best L1 minimization to continuous f
% discontinuous f is not supported yet but on the list
% f: function in chebfun
% n: degree of p
% xiini: initial interpolation points (optional)
% reps: [rep1 rep2] do rep1 steepest descent, then rep2 Newton
% reps: [rep1 rep2 reps3] reps3 points in linprog
% 
% please use [-1,1] for the moment

doLP = 0; 

if doLP
if length(f)<=n+1 & length(f.ends)==2, % if n larger than deg(f), then trivially f=p
    disp('n larger than deg(f), f=p');
    p = f; xi = []; sumnorms = []; sumshistory = [];
    return
end

ab = f.ends;
a = ab(1); b = ab(2); 
fori = f;

if ab(1)~= -1 || ab(end)~= 1;
f = chebfun(@(x)f((b-a)/2*x+(a+b)/2),'splitting','on');
end
a = -1; b = 1;

% Try chebyshev interpolation first, if successful it gives trivial solution 
xi = chebpts(n+3,[a b]); xi = xi(2:end-1); % these are chebpts(n+3,[a,b]) setminus {a,b} by AT, more stable
ptmp = chebfun.interp1(xi,f(xi),[a b]); % Theorem 14.4 and 14.5 of Powell says pn = interpolant of f at xi. 
%xii = .5*(a+b) + .5*(b-a)*cos((n+1-(0:n))*pi/(n+2));  % these are chebpts(n+3,[a,b]) setminus {a,b}
%ptmp2 = chebfun.interp1(xii,f(xii),[a b]); norm(ptmp-ptmp2) % Theorem 14.4 and 14.5 of Powell says pn = interpolant of f at xi. 

r = roots(f - ptmp);

%{
% check sign change
rmid = (r(1:end-1)+r(2:end))/2;
rmid = [a;rmid;b];
errmid = f(rmid)-ptmp(rmid); 
rsize = 0; 
badpos = []; % remove these r
for ii = 1:length(errmid)-1
    if errmid(ii)*errmid(ii+1)<0 % sign change
    rsize = rsize+1;
    else
    badpos = [badpos;ii];
    end    
end
r(badpos+1) = []; % remove bad guys 
%}

%{
r = roots(f - ptmp,'complex');
% merge complex roots into real
ix = find(abs(imag(r))>0 & abs(imag(r))<tol); 
diffr = diff(r);
ix = ix(1:2:end); 
r(ix) = real(r(ix)+r(ix+1))/2;
r(ix+1) = [];
%}

%{
dr = diff(f - ptmp); 

tol = 1000*sqrt(eps);

%{
for ixx = 1:length(ix)
if abs(dr(r(ix)))<tol
r(ix) = real(r(ix)+r(ix+1))/2;
r(ix+1) = nan;
end
end
%}
r = sort(real(r),'ascend'); 

% merge spurious-double roots 
diffr = diff(r); ix = find(diffr<tol); 
r(ix) = (r(ix)+r(ix+1))/2; r(ix+1) = [];

%{
if abs(dr(r(ix)))<tol
r(ix) = (r(ix)+r(ix+1))/2;
r(ix+1) = [];
end
%}
%}

if length(r) == length(xi), 
    disp('Chebypoints are optimal, easy case')
    p = ptmp;
    rootstmp = r;
    sumnorms = zeros(1,length(xi));
    sumshistory = sumnorms;
    return
end
    disp('Chebypoints are NOT optimal turn to LP')

% choose size for Linear Programming (larger means slower but Newton will converge surer and faster)
%LPpoints =  min(1000+50*n,10000);
LPpoints =  min(1000+50*n,5000);
%LPpoints = LPpoints*2;
%LPpoints =  min(1000+50*n);
%LPpoints =  min(200+30*n+min(10*length(f),10000));
%LPpoints = 5000;
%LPpoints =     max(3000+50*n);
%LPpoints =  min(200+20*n,5000);
%LPpoints = 400;

% Ok the problem is the nontrivial case. Do LP followed by Newton. 

if nargin>2
rep1 = reps(1);rep2 = reps(2);
if length(reps)>2, 
LPpoints = reps(3);
end
else
    rep1 = 0; rep2 = 10;  % maximum iterations  for steepest descend + Newton 
end


if nargin>2 & isempty(xiini)==0
ptmp = chebfun.interp1(xiini,f(xiini),[a b]); % Theorem 14.4 and 14.5 of Powell says pn = interpolant of f at xi. 
end

% use discrete ell1 solution for initial guess
%[ptmp,m,ptmpini] = L1discretequad(f,LPpoints,n+1); % re-mesh for O(1/n^2) convergence (instead of O(1/n))
[ptmp,m,ptmpini] = L1discretequadcheb(f,LPpoints,n+1); % re-mesh for O(1/n^2) convergence (instead of O(1/n))

% check for corrupted polynomial
tolc = 1e-8; % tolerance
if sum(ptmp>f-tolc)+sum(ptmp<f+tolc)>2+1e-2, % f=p in set of nonzero measure
    disp('possibly a corrupted polynomial')    
    %{
    %verify this 
    r1 = roots(ptmp-f-tolc); 
    r2 = roots(ptmp-f+tolc); 
    sum(diff(r1)) + sum(diff(r2))
    %}
    p = ptmp; 
    rootstmp = [];    sumnorms = [];    sumshistory = [];    
    return
end

%ptmp = L1discrete(f,LPpoints,n+1); % O(1/n)


r = roots(f-ptmpini);
sumshistory = norm(sign(f(-1)-ptmpini(-1))*intpsign(length(xi)-1,r));

r = roots(f-ptmp);
sums = sign(f(-1)-ptmp(-1))*intpsign(length(xi)-1,r);
sumshistory = [sumshistory norm(sums)];

fprintf('LPerr first second:  %9.2e,  %9.2e\n',     sumshistory)
end % doLP

if exist('ptmp','var')==0
a = -1; b = 1;
xi = .5*(a+b) + .5*(b-a)*cos((n+1-(0:n))*pi/(n+2));  % these are chebpts(n+3,[a,b]) setminus {a,b}
ptmp = chebfun.interp1(xi,f(xi),[a b]); % Theorem 14.4 and 14.5 of Powell says pn = interpolant of f at xi. 
r = roots(f-ptmp);
end

sums = sign(f(-1)-ptmp(-1))*intpsign(length(xi)-1,r);
sumshistory = norm(sums);

itwatson = 30; 
for ii = 1:itwatson % Watson update
g = sign(f(-1)-ptmp(-1))*intpsign(length(xi)-1,r);

T = chebpoly(0:n);
A = T(r,:);
dfp = diff(f-ptmp);
D = zeros(1,length(r));
for jj = 1:length(r)
    D(jj) = dfp(r(jj));
end
%D = diag(abs(D))/2; H = A'*inv(D)*A;
D = diag(2./abs(D)); H = A'*D*A;
if cond(H)>1e10, keyboard, H = H+norm(H)*eye(size(H))*1e-8; end
dc = H\(g');

dp = chebfun(dc,'coeffs'); % update poly

gam = 1; 
pnow = ptmp; 
while gam>1e-5
ptmp = pnow + gam*dp;
r = roots(f-ptmp);
sums = sign(f(-1)-ptmp(-1))*intpsign(length(xi)-1,r);
if norm(sums)<sumshistory(end) & length(r)>=n+1, break, end
    gam = gam*.8; 
    [gam norm(sums)]
end
sumshistory = [sumshistory norm(sums)]

if sumshistory(end)<1e-10, 'L1min converged', break; end
end
p = ptmp;
return



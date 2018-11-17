function F = BMCIII(f,m,n,p)
% Take a function handle and a size and return the coefficients in the
% CFF basis by evaluating it at only a fourth of the domain

% Evaluation points (assuming m odd, n even and p even > 4)
r = reshape(chebpts(m),m,1,1);
lam = reshape([pi*trigpts(n); pi],1,n+1,1);
th = reshape([pi*trigpts(p); pi],1,1,p+1);

% Array of coefficients
F = zeros(m,n,p);

% Add dependence on r, lambda and theta
f1 = @(r,lam,th)f(r,lam,th) + 0*r + 0*lam + 0*th;

%% g : evaluation at [0,1] x [-pi,0] x [0,pi]
g = feval(f1, r(floor(m/2)+1:m), lam(1:n/2+1), th(p/2+1:p+1));

%% h : evaluation at [0,1] x [0,pi] x [0,pi]
h = f1(r(floor(m/2)+1:m), lam(n/2+1:end), th(p/2+1:p+1));

%% Impose BMC-III Structure
% f(0,:,:) = constant
% Count lambda = pi only once
% Compute the mean of f evaluated at r = 0
m_zeroR = mean([mean2(g(1,:,:)), mean2(h(1,2:end,:))]);
g(1,:,:) = m_zeroR;
h(1,:,:) = m_zeroR;

% f(r,:,0) = constant
m_zeroT = mean([mean(g(:,:,1),2),mean(h(:,2:end,1),2)],2);
g(:,:,1) = repmat(m_zeroT,1,size(g,2));
h(:,:,1) = repmat(m_zeroT,1,size(g,2));

% f(r,:,pi) = constant
m_piT = mean([mean(g(:,:,end),2),mean(h(:,2:end,end),2)],2);
g(:,:,end) = repmat(m_piT,1,size(g,2));
h(:,:,end) = repmat(m_piT,1,size(g,2));

%% Flip g and h on the radial direction
flip1g = flip(g(1+mod(m,2):end,:,:), 1);
flip1h = flip(h(1+mod(m,2):end,:,:), 1);

%% Fill in F
% [0,1] x [-pi,0] x [0,pi[
F(floor(m/2)+1:m, 1:n/2+1, floor((p+1)/2)+1:p) = g(:,:,1:end-1);
% [0,1] x [0,pi[ x [0,pi[
F(floor(m/2)+1:m, n/2+1:n, floor((p+1)/2)+1:p) = h(:,1:end-1,1:end-1);
% [-1,0[ x [-pi,0] x [0,pi[
F(1:floor(m/2), 1:n/2+1, floor((p+1)/2)+1:p) = flip(flip1h(:,:,2:end),3);
% [-1,0[ x [0,pi[ x [0,pi[
F(1:floor(m/2), n/2+1:n, floor((p+1)/2)+1:p) = flip(flip1g(:,1:end-1,2:end),3);
% [0,1] x [-pi,0] x [-pi,0]
F(floor(m/2)+1:m, 1:n/2+1, 1:floor((p+1)/2)) = flip(h(:,:,2:end),3);
% [0,1] x [0,pi[ x [-pi,0]
F(floor(m/2)+1:m, n/2+1:n, 1:floor((p+1)/2)) = flip(g(:,1:end-1,2:end),3);
% [-1,0[ x [0,pi[ x [-pi,0]
F(1:floor(m/2), n/2+1:n, 1:floor((p+1)/2)) = flip1h(:,1:end-1,1:end-1);
% [-1,0[ x [-pi,0] x [-pi,0]
F(1:floor(m/2), 1:n/2+1, 1:floor((p+1)/2)) = flip1g(:,:,1:end-1);
end
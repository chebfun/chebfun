function F = BMCIII(f,m,n,p)
% Take a function handle and a size and return the coefficients in the
% CFF basis by evaluating it at only a fourth of the domain

% Evaluation points
r = reshape(chebpts(m),m,1,1);
lam = reshape([pi*trigpts(n); pi],1,n+1,1);
th = reshape([pi*trigpts(p); pi],1,1,p+1);

% Array of coefficients
F = zeros(m,n,p);

f1 = @(r,lam,th)f(r,lam,th) + 0*r + 0*lam + 0*th;

% We assume n is even

%% Chebyshev
% m even : evaluate at half and then flip it
% m odd : evaluate at floor(m/2)+1 and then flip floor(m/2)+1:m
% -> Evalute at r(floor(m/2)+1:m) then flip r(floor(m/2)+1+mod(m,2):m)

%% Fourier Theta
% Evaluate at th(floor((p+1)/2)+1:p+1) then 
% [flip(th(floor((p+1)/2)+1:p+1)) ; th(floor((p+1)/2)+1:p)]

%% g
% Evaluate at r(floor(m/2)+1:m), lam(1:n/2+1), th(floor((p+1)/2)+1:p+1)
g = f1(r(floor(m/2)+1:m), lam(1:n/2+1), th(floor((p+1)/2)+1:p+1));

%% h
% Evaluate at r(floor(m/2)+1:m), lam(n/2+1:end), th(floor((p+1)/2)+1:p+1)
h = f1(r(floor(m/2)+1:m), lam(n/2+1:end), th(floor((p+1)/2)+1:p+1));

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
F(floor(m/2)+1:m, 1:n/2+1, 1:floor((p+1)/2)) = flip(h(:,:,1+mod(p+1,2):end),3);
% [0,1] x [0,pi[ x [-pi,0]
F(floor(m/2)+1:m, n/2+1:n, 1:floor((p+1)/2)) = flip(g(:,1:end-1,1+mod(p+1,2):end),3);
% [-1,0[ x [0,pi[ x [-pi,0]
F(1:floor(m/2), n/2+1:n, 1:floor((p+1)/2)) = flip1h(:,1:end-1,1:end-mod(p+1,2));
% [-1,0[ x [-pi,0] x [-pi,0]
F(1:floor(m/2), 1:n/2+1, 1:floor((p+1)/2)) = flip1g(:,:,1:end-mod(p+1,2));
end
function F = BMCIII(f,m,n,p)
% Take a function handle and a size and return the coefficients in the
% CFF basis by evaluating it at only a fourth of the domain

% Remove issues when f is constant
f1 = @(r,lam,th) f(r,lam,th) + 0*r + 0*lam + 0*th;

% Evaluation points
r = chebpts(m);
lam = [pi*trigpts(n); pi];
th = [pi*trigpts(p); pi];

% Array of coefficients
F = zeros(m,n,p);

%% Chebyshev
% m even : evaluate at half and then flip it
% m odd : evaluate at floor(m/2)+1 and then flip floor(m/2)+1:m
% -> Evalute at r(floor(m/2)+1:m) then flip r(floor(m/2)+1+mod(m,2):m)

%% Fourier Theta
% Evaluate at th(floor((p+1)/2)+1:p+1) then 
% [flip(th(floor((p+1)/2)+1:p+1)) ; th(floor((p+1)/2)+1:p)]

%% g
% Evaluate at r(floor(m/2)+1:m), lam(1:floor(n/2)), th(floor((p+1)/2)+1:p+1)
g = reshape(f1(r(floor(m/2)+1:m), lam(1:floor((n+1)/2)), th(floor((p+1)/2)+1:p+1)),a,b,c);

%% h
% Evaluate at r(floor(m/2)+1:m), lam(floor(n/2)+1:n), th(floor((p+1)/2)+1:p+1)
h = reshape(f1(r(floor(m/2)+1:m), lam(floor((n+1)/2)+1:n+1), th(floor((p+1)/2)+1:p+1)),a,b,c);

%% Flip g and h on the radial direction
flip1g = flip(g(1+mod(m,2):end,:,:), 1);
flip1h = flip(h(1+mod(m,2):end,:,:), 1);

%% Fill in F
% [0,1] x [-pi,0] x [0,pi]
F(floor(m/2)+1:m, 1:floor((n+1)/2), floor((p+1)/2)+1:p) = g(:,:,1:end-1);
% [0,1] x [0,pi] x [0,pi]
F(floor(m/2)+1:m, floor((n+1)/2)+1:n, floor((p+1)/2)+1:p) = h(:,1:end-1,1:end-1);
% [-1,0] x [-pi,0] x [0,pi]
F(1:floor(m/2), 1:floor((n+1)/2), floor((p+1)/2)+1:p) = flip1g(:,:,1:end-1);
% [-1,0] x [0,pi] x [0,pi]
F(1:floor(m/2), floor((n+1)/2)+1:n, floor((p+1)/2)+1:p) = flip1h(:,1:end-1,1:end-1);
% [0,1] x [-pi,0] x [-pi,0]
F(floor(m/2)+1:m, 1:floor((n+1)/2), 1:floor((p+1)/2)) = flip(h,3);
% [0,1] x [0,pi] x [-pi,0]
F(floor(m/2)+1:m, 1:floor((n+1)/2), 1:floor((p+1)/2)) = flip(g,3);
% [-1,0] x [0,pi] x [-pi,0]
F(1:floor(m/2), floor((n+1)/2)+1:n, 1:floor((p+1)/2)) = flip(flip1h(:,1:end-1,1+mod(p+1,2):end),3);
% [-1,0] x [-pi,0] x [-pi,0]
F(1:floor(m/2), 1:floor((n+1)/2), 1:floor((p+1)/2)) = flip(flip1g(:,:,1+mod(p+1,2):end),3);
end
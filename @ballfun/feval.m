function H = feval(f,r,lambda,theta)
% FEVAL Evaluate a BALLFUN function
%   FEVAL(f, r, lambda, theta) is the array of values of the BALLFUN
%   function f at the grid r x lambda x theta

F = f.coeffs;
[m,n,p] = size(f);

% Get the size of the lists
Nr = length(r);
Nlam = length(lambda);
Nth = length(theta);

% Transform the lists to vectors
r = reshape(r, Nr, 1);
lambda = reshape(lambda, Nlam, 1);
theta = reshape(theta, Nth, 1);

% Fourier functions evaluated at theta
Flambda = exp(1i*lambda.*((1:n)-floor(n/2)-1));

% Fourier functions evaluated at lambda
Ftheta = exp(1i*theta.*((1:p)-floor(p/2)-1));

% Chebyshev functions evaluated at r
T = zeros(Nr,m);
T(:,1) = ones(Nr,1); T(:,2) = r;
for i = 3:m
    T(:,i) = 2*r.*T(:,i-1)-T(:,i-2);
end

G = zeros(Nr,Nlam,p);
% Evaluate f at the points r and lambda
for i = 1:p
    G(:,:,i) = T*F(:,:,i)*Flambda.';
end

% Permute G to evaluate f at theta
G = permute(G, [3,1,2]);

H = zeros(Nth,Nr,Nlam);
% Evaluate f at the points theta
for i = 1:Nlam
   H(:,:,i) = Ftheta*G(:,:,i); 
end

% Permute H to get the array of values r x lambda x theta
H = permute(H, [2,3,1]);
end

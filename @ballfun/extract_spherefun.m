function g = extract_spherefun(f, varargin)
% EXTRACT_SPHEREFUN SPHEREFUN corresponding to the value of f at radius r
%   EXTRACT_SPHEREFUN(f, r) is the SPHEREFUN function 
%   g(lambda, theta) = f(r, :, :)

F = f.coeffs;
[m,n,p] = size(f);

% Find the radius r
if nargin > 1
    r = varargin{1};
else
    r = 1;
end

if m == 1
    G = reshape(F(1,:,:),n,p);
else
    % Chebyshev functions evaluated at r
    T = zeros(1,m);
    T(1) = 1; T(2) = r;
    for i = 3:m
        T(i) = 2*r*T(i-1)-T(i-2);
    end

    % Build the array of coefficient of the spherefun function
    G = zeros(n,p);
    for i = 1:p
        G(:,i) = T*F(:,:,i);
    end
end
% Build the spherefun function; coeffs2spherefun takes the theta*lambda matrix
% of coefficients
g = spherefun.coeffs2spherefun(G.');
end

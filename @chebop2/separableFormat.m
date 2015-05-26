function [cellU, S, cellV] = separableFormat(A, xorder, yorder, domain)
%SEPARATBLEFORMAT  Compute separable expression for a linear PDO.
%
% Calculate a separable representation of a partial differential 
% operator. These representations can then be using to derive a 2D spectral
% method from 1D ideas.  The linear PDO can have variable coefficients. 
%
% This uses the tensor-train decomposition and Proposition
% 4.2 from [1].
% 
% [1] A. Townsend and S. Olver, The automatic solution of partial differential
% equations using a global spectral method, submitted, 2014. 
% 
% Author: Alex Townsend September 2014.

if ( nargin == 1 )
    N = A; 
    A = N.coeffs; 
    xorder = N.xorder; 
    yorder = N.yorder; 
    domain = N.domain; 
end

% Loop over coefficients of A. Find coefficient of highest degree. We will 
% need this to recover the variable coefficients:  
n = 10; 
for jj = 1:size(A, 1) 
    for kk = 1:size(A, 2); 
        if ( isa(A{jj,kk}, 'chebfun2') )
            [xdeg, ydeg] = length(A{jj,kk});    % get degrees 
            n = max([xdeg, ydeg, n]) + 1;       % take maximum degree we find
        end
    end
end

% Set up Chebyshev points that we need:  
x = chebpts(xorder+1); 
s = chebpts(n, domain(1:2)); 
y = chebpts(yorder+1); 
t = chebpts(n, domain(3:4)); 
[xx, ss, yy, tt] = ndgrid(x, s, y, t); 
[newx, newy] = meshgrid(s, t); 

% We need to apply Proposition 4.2 from [1]. Use the linear operator 
% T motivated by umbral calculus. That is, convert 
% 
%       d^j/dx^j-> x^j 
%       d^j/dy^j-> y^j
%         x     -> s
%         y     -> t 
% We obtain a function of 4 variables, H(x,s,y,t). See [1]. 
H = @(x,s,y,t) 0*x;
for jj = 1:size(A, 1)
    for kk = 1:size(A, 2)
        if ( isa(A{jj,kk}, 'double') && ~(A{jj,kk} == 0) )
            H = @(x,s,y,t) H(x,s,y,t) + A{jj,kk}*x.^(kk - 1).*y.^(jj - 1);
        elseif ( isa(A{jj,kk}, 'chebfun2') )
            v = zeros(1, n, 1, n);
            v(1,:,1,:) = feval(A{jj,kk}, newx, newy).';
            out = repmat(v, [xorder + 1, 1, yorder + 1, 1]);
            H = @(x,s,y,t) H(x,s,y,t) + out.*x.^(kk-1).*y.^(jj-1);
        end
    end
end
H = H(xx, ss, yy, tt);

% Using tensor-train ideas. Calculate the splitting rank of the function
% H(x,s,y,t): 
A = reshape(H, n*(xorder+1), n*(yorder+1));
[U, S, V] = svd(A); 
rk = find(abs(diag(S)/S(1,1)) > 1000*eps, 1, 'last' );  % splitting rank

% Restrict to singular vectors of interest.  
S = S(1:rk, 1:rk);
U = U(:, 1:rk); 
V = V(:, 1:rk); 

% We have the splitting rank of the PDO, now we want the corresponding 
% separable representation. The following is tricky to get right... 
cellU = cell(yorder+1, rk);
cellV = cell(xorder+1, rk); 
c1 = cell(xorder, 1); 
c2 = cell(yorder, 1); 
% Matrices to convert ChebT -> monomials: 
converty = fliplr(poly(chebpoly(0:yorder))); 
convertx = fliplr(poly(chebpoly(0:xorder)));

% Figure out the separable representation: 
for jj = 1:rk 

    % This is giving us the 1D ODEs that go on the right in the generalized
    % Sylvester matrix equation: 
    f1 = chebfun2.vals2coeffs( reshape(U(:,jj), xorder+1, n) );  % @(x, s)
    f1 = f1.' * convertx;
    for kk = 1:xorder+1
        c1{kk} = chebfun(f1(:,kk), domain(1:2), 'coeffs');
        cellV(kk,jj) = c1(kk);
    end

    % This is giving us the 1D ODEs that go on the left in the generalized
    % Sylvester matrix equation: 
    f2 = chebfun2.vals2coeffs( reshape(conj(V(:,jj)), yorder+1, n) );  % @(y, t) 
    f2 = f2.' * converty;
    for kk = 1:yorder+1
        c2{kk} = chebfun(f2(:,kk), domain(3:4), 'coeffs');
        cellU(kk,jj) = c2(kk);
    end  
end

end
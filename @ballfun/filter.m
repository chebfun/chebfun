function g = filter( f )
% FILTER Filter of a nonsmooth BALLFUN function
%   FILTER(f) removes the discontinuities of the BALLFUN function f

% Copyright 2018 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Test with f(x,y,z) = cos(x*y):
% f=ballfun(@(r,lam,th)cos(r^2*sin(th)^2*cos(lam)*sin(lam)),[20,20,20]);

[m,n,p] = size(f);

% Coefficients of f
F = f.coeffs;
% Permute F in r x th x lam to transform the coefficients to values
F = permute(F, [1,3,2]);

% Transform F into Vals_R x Vals_Th x Fourier_Lam
for j = 1:n
    F(:,:,j) = chebtech2.coeffs2vals( F(:,:,j) );
    F(:,:,j) = trigtech.coeffs2vals( F(:,:,j).' ).';
end

% Tensor mask for X in Vals_Th x Fourier_Lam x Vals_R
Tmask = zeros(p,n,m);

% Get the mask of the disk
mask = DiskFilter(m,n);

% Fill in Tmask(:,:,i) with the filter of the sphere (take into account the
% radius associated with i)
for i = 1:m
   gamma = sum(mask(i,1:floor(n/2)+1));
   Mr = maskR(p,n,gamma);
   Tmask(:,:,i) = Mr;
end

% Permute Tmask to the system Vals_R x Vals_Th x Fourier_Lam 
Tmask = permute(Tmask,[3,1,2]);

% Apply the filter to F
G = F.*Tmask;

% Transform it to an array of coeffs
for j = 1:n
    G(:,:,j) = chebtech2.vals2coeffs( G(:,:,j) );
    G(:,:,j) = trigtech.vals2coeffs( G(:,:,j).' ).';
end

% Permute back
G = permute(G, [1,3,2]);

% Return the ballfun function associated with G
g = ballfun(G);

end

% Build the filter on the disk
function mask = DiskFilter(m,n)
mask = zeros( m, n ); 
mask( :, floor(n/2) + 1 ) = 1; 
k = 1:floor(n/2);
j = 1:floor(m/2);
[kk,jj] = meshgrid( k, j );
% Fill in one triangle
mask( j, k ) = ( 0.9*jj < kk ); 
% Do permutations to fill in the matrix
mask( j, floor(n/2)+2:n ) = fliplr( mask( j,k(2-mod(n,2):floor(n/2)) ) );
mask( m-floor(m/2)+1:m, :) = flipud( mask( j(1:floor(m/2)), :) ); 
end

% Build the mask th*lam for a fixed r (half of the width of the
% parallelogram = gamma
function Mr = maskR(P,N,gamma) 
mask = zeros( P, N );
n=N; p=floor(P/2);
mask( floor(p/2)+1, floor(n/2) + 1 ) = 1; 
k = 1:floor(n/2)+1;
j = 1:floor(p/2);

% Compute the equation of the line th = g(lam) = a*lam + b 
% slope
a = (floor(p/2)-1)/(floor(n/2)+1-(floor(n/2)-gamma+2));
% gamma corresponds to the half-width of the parallelogram
b = 1-(floor(n/2)-gamma+2)*a;
[kk,jj] = meshgrid( k, j );

% Fill in one triangle
if gamma ~= 1
    % change it to ( jj <= ceil(a*kk+b)) to add more tolerance
    mask( j, k ) = ( jj <= a*kk+b);
else
   mask(j,floor(n/2)+1) = 1; 
end

% Do permutations to fill in the matrix
mask( j, floor(n/2)+2:n ) = fliplr( mask( j,k(2-mod(n,2):floor(n/2)) ) );
mask( p-floor(p/2)+1-mod(p,2):p, :) = mask(1:floor(p/2)+mod(p,2), :);
mask( 1:floor(p/2), :) = flipud( mask( p-floor(p/2)+1:p, :) ); 
mask( P-floor(P/2)+1:P, :) = mask(1:floor(P/2), :);

% Fill in the central coefficient
mask( floor(P/2)+1, floor(N/2) + 1 ) = 1; 

% Return the mask
Mr = mask;

end

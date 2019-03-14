function X = coeffs2vals( X )
%COEFFS2VALS Convert an array of Chebyshev--Fourier--Fourier of
%            coefficients to an array of values.
%
%   COEFFS2VALS( CFS ) computes the values on a tensor-product doubled-up
%   spherical grid of a function that is represented by a
%   Chebyshev-Fourier-Fourier expansion with expansion coefficients CFS.
%
% See also CHEBFUN.COEFFS2VALS, VALS2COEFFS.

% Copyright 2019 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if isempty( X )
    return
end

[m, n, p] = size( X );

% Old slow approach:
% for k = 1:p
%     X(:,:,k) = chebtech2.coeffs2vals( X(:,:,k) );
%     X(:,:,k) = trigtech.coeffs2vals( X(:,:,k).' ).';
% end
% for j = 1:n
%     vj = reshape( X(:,j,:), m, p );
%     vj = trigtech.coeffs2vals( vj.' ).';
%     X(:,j,:) = reshape( vj, m, 1, p );
% end

% Faster approach, but code less readable:
if ( m > 1 ) 
    X(2:m-1, :, :) = X(2:m-1, :, :)/2; 
    X = fft( vertcat(X, X(m-1:-1:2,:,:)), [], 1 );
    X = X(m:-1:1, :, :);
end

scl_p = (n*p)*even_odd_fix( p );
scl_n = even_odd_fix( n );
Enp = reshape(scl_n.'*scl_p, [1 n p]);

X = ifft(ifft(ifftshift(ifftshift(X.*repmat( Enp, m, 1, 1),2),3), [], 2),[],3);

end

function scl = even_odd_fix( n )

if ( mod(n, 2) ) 
    scl = (-1).^(-(n-1)/2:(n-1)/2);
else
    scl = (-1).^((-n/2):(n/2-1));
end

end
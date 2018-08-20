function X = vals2coeffs( X )
% VALS2COEFFS Convert a Chebyshev--Fourier--Fourier array of coefficients 
% to an array of values
%   VALS2COEFFS(VALS) is the array of values corresponding to VALS

% Copyright 2018 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

[m, n, p] = size( X );

% Old, slow, approach: 
% for k = 1:p
%     X(:,:,k) = chebtech2.vals2coeffs( X(:,:,k) );
%     X(:,:,k) = trigtech.vals2coeffs( X(:,:,k).' ).';
% end
% for j = 1:n
%     vj = reshape( X(:,j,:), m, p );
%     vj = trigtech.vals2coeffs( vj.' ).';
%     X(:,j,:) = reshape( vj, m, 1, p );
% end

% Faster approach, but less readable code:
X = ifft( vertcat( X(m:-1:2,:,:),X ), 2*(m-1), 1);
X = X(1:m, :, :);
X(1,:,:) = X(1,:,:)/2;
X(m,:,:) = X(m,:,:)/2;

X = fftshift(fftshift(fft(fft(X, [], 2),[],3), 2), 3);

scl_p = (2/n/p)*even_odd_fix( p );
scl_n = even_odd_fix( n );
Enp = reshape(scl_n.'*scl_p, [1 n p]);
X = X.*repmat( Enp, m, 1, 1);

end

function scl = even_odd_fix( n )

if ( mod(n, 2) ) 
    scl = (-1).^(-(n-1)/2:(n-1)/2);
else
    scl = (-1).^((-n/2):(n/2-1));
end

end
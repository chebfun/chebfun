function X = vals2coeffs( X )
%VALS2COEFFS   BALLFUN values to coefficients.
%
% C = VALS2COEFFS( V ) converts values in (r, lambda, theta)
% in spherical coordinates at a chebpts x trigpts x trigpts grid to
% Chebyshev x Fourier x Fourier coefficients.
%
% This command is mainly for internal use.
%
% See also CHEBFUN.VALS2COEFFS, COEFFS2VALS.

% Copyright 2019 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if isempty( X )
    return
end

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

if ( m > 1 )
    X = ifft( vertcat( X(m:-1:2,:,:),X ), 2*(m-1), 1);
    X = 2*X(1:m, :, :);
    X(1,:,:) = X(1,:,:)/2;
    X(m,:,:) = X(m,:,:)/2;
end

X = fftshift(fftshift(fft(fft(X, [], 2),[],3), 2), 3);

scl_p = (1/n/p)*even_odd_fix( p );
scl_n = even_odd_fix( n );
Enp = reshape(scl_n.'*scl_p, [1 n p]);
X = X.*repmat( Enp, m, 1, 1);

end

function scl = even_odd_fix( n )
% Even/odd scaling:
if ( mod(n, 2) ) 
    scl = (-1).^(-(n-1)/2:(n-1)/2);
else
    scl = (-1).^((-n/2):(n/2-1));
end

end
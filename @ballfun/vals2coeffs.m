function X = vals2coeffs( X )
%VALS2COEFFS Convert an array of values to an array of coefficients in a 
%            Chebyshev--Fourier--Fourier expansion.
%
%   C = VALS2COEFFS( F ) computes the Chebyshev-Fourier-Fourier expansion
%   coefficients for the array of values F.  F is assumed to be sampled on
%   a doubled-up spherical grid of m-n-p points, where m corresponds to the
%   radial sampling at Chebyshev points from [-1,1], n corresponds to
%   azimuthal sampling at equally spaced points from (-pi,pi], and p
%   corresponds to polar sampling at equally spaced points from (-pi,pi].
% 
% See also CHEBFUN.COEFFS2VALS, COEFFS2VALS.

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
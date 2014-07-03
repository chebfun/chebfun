function D = diffmat(N, m)
%DIFFMAT   Fourier differentiation matrix.
%   D = DIFFMAT(N) is the matrix that maps function values at N equally-spaced
%   points to values of the derivative of the interpolating trigonometric
%   polynomial at those points.
%
%   D = DIFFMAT(N, K) is the same, but for the Kth derivative.
%
%   The matrices are computed using the formulae given in Spectral methods in
%   Matlab [1].

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% References:
%  [1] L.N. Trefethen, "Spectral Methods in Matlab", SIAM, Philadelphia, 2000. 

% Grid point spacing h.
h = 2*pi/N; 

% No differentiation.
if ( m == 0 )
    D = eye(n);

% First-order Fourier differentiation matrix.
elseif ( m == 1 )
    
    if ( mod(N, 2) ) % N is odd.
        column = [0 .5*csc((1:N-1)*h/2)]';
    else % N is even.
        column = [0 .5*cot((1:N-1)*h/2)]';
    end
    column(2:2:end) = -column(2:2:end);
    row = column([1 N:-1:2]);
    D = toeplitz(column, row);

% Second-order Fourier differentiation matrix.
elseif ( m == 2 )
    
    if ( mod(N, 2) ) % N is odd.
        tmp = csc((1:N-1)*h/2).*cot((1:N-1)*h/2);
        column = [pi^2/3/h^2-1/12 .5*tmp].';
    else % N is even.
        column = [pi^2/3/h^2+1/6 .5*csc((1:N-1)*h/2).^2].';
    end
    column(1:2:end) = -column(1:2:end);
    D = toeplitz(column);
    
% Higher-orders Fourier differentiation matrices.
else
    
    % [TODO]: Improve efficiency of this code for higher derivatives.
    if ( mod(N, 2) ) % N is odd.
        column = (1i*[0:(N-1)/2 -(N-1)/2:-1]').^m;
    else % N is even.
        column = (1i*[0:N/2-1 0 -N/2+1:-1]').^m;
    end
    D = real(ifft(bsxfun(@times, column, fft(eye(N)))));
    
end

end
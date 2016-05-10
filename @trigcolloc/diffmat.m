function D = diffmat(N, m)
%DIFFMAT   Trigonometric Fourier differentiation matrix.
%   D = DIFFMAT(N) is the matrix that maps function values at N equally-spaced
%   points in [-1, 1) to values of the derivative of the interpolating
%   trigonometric polynomial at those points.
%
%   D = DIFFMAT(N, K) is the same, but for the Kth derivative.
%
%   The matrices are computed using the formulae given in [1].
%
% References:
%   [1] L.N. Trefethen, "Spectral Methods in MATLAB". SIAM, Philadelphia, 2000.
%   [2] J.A.C. Weideman and S.C. Reddy, A MATLAB differentiation matrix suite, 
%       ACM TOMS, Vol. 26, No. 4, Page 465--519, 2000. 

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% First deal with some trivial cases.
if ( N == 0 )
    D = [];
    return
end
if ( N == 1 )
    D = 0;
    return
end

% Default differentiation order is 1:
if ( nargin < 2 )
    m = 1;
end

% Grid point spacing h. (Note that wee start off on [-pi,pi) and scale to [-1,1)
% as the final step.
h = 2*pi/N;

% No differentiation: identity matrix.
if ( m == 0 )
    D = eye(N);
    return
    
    % First-order trigonometric Fourier differentiation matrix.
elseif ( m == 1 )
    
    if ( mod(N, 2) ) % N is odd.
        column = [0, .5*csc((1:N-1)*h/2)]';
    else % N is even.
        column = [0, .5*cot((1:N-1)*h/2)]';
    end
    column(2:2:end) = -column(2:2:end);
    row = column([1 N:-1:2]);
    D = toeplitz(column, row);
    
    % Second-order trigonometric Fourier differentiation matrix.
elseif ( m == 2 )
    
    if ( mod(N, 2) ) % N is odd.
        tmp = csc((1:N-1)*h/2).*cot((1:N-1)*h/2);
        column = [pi^2/3/h^2 - 1/12, .5*tmp].';
    else % N is even.
        column = [pi^2/3/h^2 + 1/6, .5*csc((1:N-1)*h/2).^2].';
    end
    column(1:2:end) = -column(1:2:end);
    D = toeplitz(column);
    
    % Third-order trigonometric Fourier differentiation matrix.
elseif ( m == 3 )
    
    if ( mod(N, 2) ) % N is odd.
        cscc = csc((1:N-1)*h/2);
        cott = cot((1:N-1)*h/2);
        column = [0, 3/8*cscc.*cott.^2 + 3/8*cscc.^3 - pi^2/2/h^2*cscc];
    else % N is even.
        tmp = csc((1:N-1)*h/2).^2.*cot((1:N-1)*h/2);
        column = [0, 3/4*tmp - pi^2/2/h^2*cot((1:N-1)*h/2)].';
    end
    column(2:2:end) = -column(2:2:end);
    row = column([1 N:-1:2]);
    D = toeplitz(column, row);
    
    % Fourth-order trigonometric Fourier differentiation matrix.
elseif ( m == 4 )
    
    cscc = csc((1:N-1)*h/2);
    cott = cot((1:N-1)*h/2);
    if ( mod(N, 2) ) % N is odd.
        column = [- pi^4/5/h^4 + pi^2/6/h^2 - 7/240, ...
            5/4*cscc.^3.*cott + 1/4*cscc.*cott.^3 - ...
            (pi^2/h^2)*cscc.*cott].';
    else % N is even.
        column = [- pi^4/5/h^4 - pi^2/3/h^2 + 1/30, ...
            cscc.^2.*cott.^2 + .5*cscc.^4 - ...
            (pi^2/h^2)*cscc.^2].';
    end
    column(1:2:end) = -column(1:2:end);
    D = toeplitz(column);
    
    % Higher-orders trigonometric Fourier differentiation matrices.
else
    
    % [TODO]: Improve efficiency of this code for higher derivatives.
    if ( mod(N, 2) ) % N is odd.
        column = (1i*[0:(N-1)/2, -(N-1)/2:-1]').^m;
    else % N is even.
        column = (1i*[0:(N/2-1), 0, (-N/2+1):-1]').^m;
    end
    D = real(ifft(bsxfun(@times, column, fft(eye(N)))));
    
end

% Scale from [-pi, pi) to [-1, 1):
D = pi^m*D;

end

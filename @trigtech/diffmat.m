function D = diffmat(n, p)
%DIFFMAT  Trigonometric differentiation matrix.
%   D = DIFFMAT(N) is the matrix that maps function values at N equally-spaced 
%   points in [0 2*pi) to values of the derivative of the Fourier interpolant 
%   at those points.
%
%   D = DIFFMAT(N, K) is the same, but for the Kth derivative.
%
% See also DIFF.

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% References for the rectangular differentiation matrices:
%
% J.A.C. Weideman and S.C. Reddy, A MATLAB differentiation matrix suite, ACM
% Transcations on Mathematical Software, Vol. 26, No. 4, Page 465--519, 2000. 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% grid spacing:
h = 2*pi/n;

% sign of entries in a column:
sgn = (-1).^(1:n-1).';

% indices:
n1 = floor((n-1)/2);
n2 = ceil((n-1)/2);

% grid points on [-1 0):
v = (1:n2)'*h/2;

% For different p:
if ( p == 0 )           % Zeroth-derivative
    D = eye(n);
else
    switch p
        case 1          % 1st-derivative
            
            % Forming the first column by 'flipping trick':
            if ( rem(n, 2) )           % n is odd
                tmp = csc(v);
                col = [0; (pi/2)*sgn.*[tmp; flipud(tmp(1:n1))]];
            else                      % n is even
                tmp = cot(v);
                col = [0; (pi/2)*sgn.*[tmp; -flipud(tmp(1:n1))]];
            end
            
            % Form the first row:
            row = -col;
            
        case 2         % 2nd-derivative
            
            % Form the first column by 'flipping trick':
            if ( rem(n, 2) )           % n is odd
                tmp = csc(v).*cot(v);
                col = pi^2*[-pi^2/(3*h^2)+1/12; ...
                    -0.5*sgn.*[tmp; -flipud(tmp(1:n1))]];
            else                      % n is even
                tmp = csc(v).^2;
                col = pi^2*[-pi^2/(3*h^2)-1/6; ...
                    -0.5*sgn.*[tmp; flipud(tmp(1:n1))]];
            end
            
            % Form the first row:
            row = col;
            
        otherwise      % pth-derivative for p >= 3
            
            % Form the first column using FFT:
            n3 = (-n/2)*rem(p+1,2)*ones(rem(n+1,2));
            waveNumber = 1i*[(0:n1) n3 (-n1:-1)];
            col = pi^p*real(ifft(waveNumber.^p.*fft(eye(1,n))));
            
            % Form the first row:
            if ( rem(p,2) )      % p is odd
                col = [0 col(2:n)]';
                row = -col;
            else                 % p is even
                row = col;
            end
    end
    
    % Form the differentiation matrix which is toeplitz:
    D = toeplitz(col, row);
    
end

end

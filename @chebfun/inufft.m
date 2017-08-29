function [c, p] = inufft( f, x, omega, type)
%CHEBFUN.INUFFT   Inverse nonuniform fast Fourier transform
%   [C, P]= CHEBFUN.INUFFT(F) is the same as ifft( F ). F must be a column
%   vector. C = P(F) is a planned version of the fast transform.
%
%   C = CHEBFUN.INUFFT(F, X) is an inverse nonuniform fast Fourier transform
%   of type 2, which computes C = A\F, where
%        A_{jk} = exp(-2*pi*1i*X(j)*k), 0<=j,K<=N-1.
%   F and X must be column vectors of the same length.
%
%   C = CHEBFUN.INUFFT(F, X, 2) is the same as CHEBFUN.INUFFT(F, X).
%
%   C = CHEBFUN.INUFFT(F, OMEGA, 1) is an inverse nonuniform fast Fourier
%   transform of type 1, which computes C=A\F, where
%        A_{jk} = exp(-2*pi*1i*j/N*OMEGA(k)), 0<=j,K<=N-1.
%
%   C = CHEBFUN.INUFFT(F, X, 1, TOL), C = CHEBFUN.INUFFT(F, OMEGA, 2, TOL),
%   and C = CHEBFUN.INUFFT(F, X, OMEGA, TOL) are the same as above but with
%   a tolerance of TOL. By default, TOL = eps.
%
%   The algorithm in this MATLAB script is based on the paper:
%       [1] D. Ruiz--Antoln and A. Townsend, "A nonuniform fast Fourier
%       transform based on low rank approximation", submitted, 2017.
%   This paper relates the NUFFT to a FFT by low rank approximation.
%   A faster MATLAB implementation outside of the Chebfun system
%   is available from the author. Please email: townsend@cornell.edu.
%
% See also chebfun.nufft and chebfun.ndct.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( nargin == 1 )
    p = @(fvals) ifft(fvals);
    c = p(f);
elseif ( nargin == 2 )
    % default to type 2 nufft
    [c, p] = inufft2( f, x, eps );
elseif ( nargin == 3 )
    type = omega;
    if ( numel(type) == 1 )
        if ( type == 1 )
            [c, p] = inufft1( f, x, eps);
        elseif ( type == 2 )
            [c, p] = inufft2( f, x, eps);
        elseif ( type<1 && type>0 && numel(f)>1 )
            tol = type;
            [c, p] = inufft1( f, x, tol);
        elseif ( numel(f) == 1 )
            error('CHEBFUN:CHEBFUN:inufft:three', ...
                'Type 3 NUIFFT has not been implemented');
        else
            error('CHEBFUN:CHEBFUN:inufft:type', ...
                'Unrecognised NUFFT type.');
        end
    elseif ( numel(type) == size(f,1) )
        % NUFFT-III:
        error('CHEBFUN:CHEBFUN:inufft:three', ...
            'Type 3 NUIFFT has not been implemented')
    else
        error('CHEBFUN:CHEBFUN:inufft:syntax',...
            'Unrecognised number of arguments to NUFFT.')
    end
elseif ( (nargin == 4) &&  (type == 3) )
        error('CHEBFUN:CHEBFUN:inufft:three', ...
            'Type 3 NUIFFT has not been implemented')
end
end

function [c, p] = inufft1(f, omega, tol)
%NUIFFT1  Compute the nonuniform IFFT of type 1.
% This is done by noting that 
%    inv(A) = A^*inv(AA^*)
%  for any matrix A.  Therefore, Ax = b can be solved in two steps: 
%   1)  Solve (A*A^*)x = b for x,
%   2)  Compute A^*x 
% When A = tilde(F)_1, then (A*A^*) is a Toeplitz matrix. 

N = size(omega,1);
[~, p]  = chebfun.nufft(f, omega, 1);
[~, pt] = chebfun.nufft(f, omega/N, 2);
pt = @(f) conj(pt(conj(f)));
col = p(ones(N,1));
row = conj(col);
row(1) = col(1); 

Afun = @(f) fastToeplitz(col, row, f);
[c, ~] = pcg(Afun, f, 100*tol, 50);
c = pt(c); 

% Planning is unavailable because of the pcg() command. 
p = [];
end

function [c, p] = inufft2(f, x, tol)
% NUIFFT2  Compute the nonuniform IFFT of type 2.
% We do this by solving the normal equations (F'*F)*f=F'*c, where F is the
% NUDFT2 matrix and using the conjugate gradient method. 

N = size(x, 1);

% tilde{F}_2^*tilde{F}_2 is a Toeplitz matrix, calculate this: 
[~, pct] = chebfun.nufft(conj(f), N*x, 1);
row = pct(ones(N,1));
pct = @(f) conj(pct(conj(f)));
col = pct(ones(N,1)); 
row(1) = col(1);

% Conjugate gradient method on normal equations:
AFun = @(c) fastToeplitz(col, row, c);
[c, ~] = pcg(AFun, pct(f), 100*tol, 50);
% Plan: 
p = @(f) pcg(AFun, pct(f), 100*tol, 50);
end

function b = fastToeplitz(col, row, x)
% Compute b = Ax, where A is a Toeplitz matrix with its first
% column given by the vector c and first row given by r.
% Note r(1) should be equal to c(1).

m = size(x, 1);

% A toeplitz matrix can be embedded into a Circulant matrix:
%  T11 = T;
%  T12 = toeplitz([0;r(end:-1:2)],[0;c(end:-1:2)]);
%  C = [T11 T12 ;
%       T12 T11];      % circulant
%  Use this fact to compute the matrix-vector multiply using FFT.
b = fastCirculant( [col ; 0 ; row(end:-1:2)], [x ; zeros(m,1)] );
b = b(1:m);

% In the 1st line, notice [v ; zeros(n,1)], which is zeroing out the
% contributions from the last n columns of C.
%
% In the 2nd line, the contributions from rows n+1:2n are removed as we do not
% need them.
%
% Therefore, we are computing a matrix-vector product only with the principal
% nxn block (the original Toeplitz matrix).

end

function b = fastCirculant(col, x)
% Compute b = Ax, where A is a circulant matrix with its first
% column given by the vector c.
d = fft( col );           % eigenvalues of A
b = ifft( d .* fft(x) );  % FFT diagonalizes A
end
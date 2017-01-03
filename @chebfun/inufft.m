function [f, p] = inufft( c, x, omega, type)
% CHEBFUN.INUFFT   Inverse nonuniform fast Fourier transform
%
% [F, P]= CHEBFUN.INUFFT( C ) is the same as ifft( C ). C must be a column
% vector. F = P(C) is a planned version of the fast transform.
%
% F = CHEBFUN.INUFFT( C, X ) is an inverse nonuniform fast Fourier transform
% of type 2, which computes A\C, where
%
%        A_{jk} = exp(-2*pi*1i*X(j)*k/N), 0<=j,K<=N-1.
%
% C and X must be column vectors of the same length.
%
% F = CHEBFUN.INUFFT( C, X, 2 ) is the same as CHEBFUN.INUFFT( C, X )
%
% F = CHEBFUN.INUFFT( C, OMEGA, 1 ) is an inverse nonuniform fast Fourier
% transform of type 1, which computes A\C, where
%
%        A_{jk} = exp(-2*pi*1i*j/N*OMEGA(k)), 0<=j,K<=N-1.
%
% F = CHEBFUN.INUFFT( C, X, 1, TOL), F = CHEBFUN.INUFFT(C, OMEGA, 2, TOL),
% and F = CHEBFUN.INUFFT( C, X, OMEGA, TOL ) are the same as above but with
% a tolerance of TOL. By default, TOL = eps.
%
% See also chebfun.nufft and chebfun.ndct.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

%%% DEVELOPER'S NOTE %%%
% The algorithm in this MATLAB script is based on the paper:
%
% [1] D. Ruiz--Antoln and A. Townsend, "A nonuniform fast Fourier transform
% based on low rank approximation", in preparation, 2016.
%
% This paper related the NUFFT to a FFT by low rank approximation.
%
% A faster MATLAB implementation outside of the Chebfun system
% is available from the author. Please email: townsend@cornell.edu.

if ( nargin == 1 )
    p = @(c) ifft(c);
    f = p(c);
elseif ( nargin == 2 )
    % default to type 2 nufft
    [f, p] = nuifft2( c, x, eps );
elseif ( nargin == 3 )
    type = omega;
    if ( numel(type) == 1 )
        if ( type == 1 )
            [f, p] = nuifft1( c, x, eps);
        elseif ( type == 2 )
            [f, p] = nuifft2( c, x, eps);
        elseif ( type<1 && type>0 && numel(c)>1 )
            tol = type;
            [f, p] = nuifft1( c, x, tol);
        elseif ( numel(c) == 1 )
            error('CHEBFUN::NUIFFT::3','Type 3 NUIFFT has not been implemented')
        else
            error('CHEBFUN::NUFFT::TYPE','Unrecognised NUFFT type.');
        end
    elseif ( numel(type) == size(c,1) )
        % NUFFT-III:
        error('CHEBFUN::NUIFFT::3','Type 3 NUIFFT has not been implemented')
    else
        error('CHEBFUN::NUFFT::SYNTAX','Unrecognised number of arguments to NUFFT.')
    end
elseif ( nargin == 4 )
    if ( type == 3 )
        error('CHEBFUN::NUIFFT::3','Type 3 NUIFFT has not been implemented')
    end
end
end

function [f,p] = nuifft1( c, omega, tol )
% NUIFFT1  Compute the nonuniform IFFT of type 1.

% This is done by noting that 
% 
%    inv(A) = A^*inv(AA^*)
% 
%  for any matrix A.  Therefore, Ax = b can be solved in two steps: 
% 
%   1)  Solve (A*A^*)x = b for x,
%   2)  Compute   A^*x 
% 
% When A = tilde(F)_1, then (A*A^*) is a Toeplitz matrix. 

N = size(omega,1);
[~, p] = chebfun.nufft( c, omega, 1);
[~, p_trans] = chebfun.nufft(c, omega/N, 2);
p_trans = @(c) conj( p_trans( conj( c ) ) );
toeplitz_col = p( ones(N,1) );
toeplitz_row = conj( p( ones(N,1) ) );
toeplitz_row(1) = toeplitz_col(1); 

toeplitz_matvec = @(c) fastToeplitz( toeplitz_col, toeplitz_row, c );
[f, ~] = pcg(toeplitz_matvec, c, 100*tol, 50 );
f = p_trans( f ); 

% Planning is unavailable because of the pcg() command. 
p = [];
end

function [f,p] = nuifft2( c, x, tol )
% NUIFFT2  Compute the nonuniform IFFT of type 2.
% We do this by solving the normal equations (F'*F)*f=F'*c, where F is the
% NUDFT2 matrix and using the conjugate gradient method. 

N = size(x,1);

% tilde{F}_2^*tilde{F}_2 is a Toeplitz matrix, calculate this: 
[~, p_ctrans] = chebfun.nufft(conj(c), N*x, 1);
toeplitz_row = p_ctrans( ones(N,1) );
p_ctrans = @(y) conj( p_ctrans( conj( y ) ) );
toeplitz_col = p_ctrans( ones(N,1) ); 
toeplitz_row(1) = toeplitz_col(1);

% Conjugate gradient method on normal equations:
toeplitz_matvec = @(c) fastToeplitz(toeplitz_col, toeplitz_row, c);
[f, ~] = pcg(toeplitz_matvec, p_ctrans(c), 100*tol, 50 );
% Plan: 
p = @(c) pcg(toeplitz_matvec, p_ctrans(c), 100*tol, 50 );
end

function b = fastToeplitz( c, r, x )
% Compute b = Ax, where A is a Toeplitz matrix with its first
% column given by the vector c and first row given by r.
% Note r(1) should be equal to c(1).

m = size( x, 1 );

% A toeplitz matrix can be embedded into a Circulant matrix:
%
%  T11 = T;
%  T12 = toeplitz([0;r(end:-1:2)],[0;c(end:-1:2)]);
%  C = [T11 T12 ;
%       T12 T11];      % circulant
%
%  Use this fact to compute the matrix-vector multiply using FFT.
b = fastCirculant( [c ; 0 ; r(end:-1:2)], [x ; zeros(m,1) ] );
b = b(1:m);

% In the 1st line, notice [v ; zeros(n,1) ], which is zeroing out the
% contributions from the last n columns of C.
%
% In the 2nd line, the contributions from rows n+1:2n are removed as we do not
% need them.
%
% Therefore, we are computing a matrix-vector product only with the principal
% nxn block (the original Toeplitz matrix).

end

function b = fastCirculant( c, x )
% Compute b = Ax, where A is a circulant matrix with its first
% column given by the vector c.

d = fft( c );             % eigenvalues of A
b = ifft( d.*fft( x ) );  % FFT diagonalizes A

end
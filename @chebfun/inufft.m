function [c, p] = inufft( f, x, type, tol)
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
%    [1] D. Ruiz-Antoln and A. Townsend, "A nonuniform fast Fourier
%    transform based on low rank approximation", SISC, 40 (2018), A529-A547.
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
    if ( numel(type) == 1 )
        if ( type == 1 )
            [c, p] = inufft1( f, x, eps);
        elseif ( type == 2 )
            [c, p] = inufft2( f, x, eps);            
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
            'Unrecognised input arguments to INUFFT.')
    end
elseif ( nargin == 4 )
    if ( numel(type) == 1 && numel(tol) == 1 )        
        if ( type == 1 )
            [c, p] = inufft1( f, x, tol);
        elseif ( type == 2 )
            [c, p] = inufft2( f, x, tol);
        else
            error('CHEBFUN:CHEBFUN:inufft:type', ...
                'Unrecognised INUFFT type.');
        end
    elseif ( numel(type) == size(f,1) )
        % NUFFT-III:
        error('CHEBFUN:CHEBFUN:inufft:three', ...
            'Type 3 NUIFFT has not been implemented')
    elseif ( numel(tol) > 1 )
        error('CHEBFUN:CHEBFUN:inufft:tol',...
            'Tolerance parameter to INUFFT must be a scalar.')
    else        
        error('CHEBFUN:CHEBFUN:inufft:syntax',...
            'Unrecognised input arguments to INUFFT.')
    end    
else
    error('CHEBFUN:CHEBFUN:inufft:syntax',...
        'Unrecognised input arguments to INUFFT.')
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
col  = chebfun.nufft(ones(N,1), omega, 1);
row = conj(col);
row(1) = col(1); 

[~, pt] = chebfun.nufft(f, omega/N, 2);
pt = @(f) conj(pt(conj(f)));

Afun = @(f) fastToeplitz(col, row, f);

    % Plan: 
    function c = plan1(f)
        [c,flag,relres] = pcg(Afun, f, 100*tol, 50);
        c = pt(c); 

        % Report a warning if the CG failed to converge
        if (flag ~= 0)
            warning('CHEBFUN:CHEBFUN:inufft1:tolerance',...
                    ['Conjugate gradient method failed to converge to a ' ...
                     'tolerance of %1.2e in %d iterations.  ' ...
                     'Tolerance reached was %1.2e.'],100*tol,50,relres);
        end    
    end

c = plan1(f);
        
p = @plan1;

end

function [c, p] = inufft2(f, x, tol)
% NUIFFT2  Compute the nonuniform IFFT of type 2.
% We do this by solving the normal equations (F'*F)*c=F'*f, where F is the
% NUDFT2 matrix and using the conjugate gradient method.

N = size(x, 1);

% tilde{F}_2^*tilde{F}_2 is a Toeplitz matrix, so we only need the first
% column, i.e. tilde{F}_2^*tilde{F}_2 e_1. Note that since tilde{F}_2 e_1 =
% ones(N,1), we really only need to compute tilde{F}_2^*ones(N,1), which 
% can be done with the NUFFT 1.
[row, pct] = chebfun.nufft(ones(N,1), N*x, 1);
col = conj(row); 
row(1) = col(1);

% Set up function for computing F'*f, which is done by applying the 
% NUFFT-I to the conjugate of f then taking the conjugate of the result.
pct = @(f) conj(pct(conj(f)));

% Conjugate gradient method on normal equations:
AFun = @(c) fastToeplitz(col, row, c);

    % Plan: 
    function c = plan2(f)
        [c,flag,relres] = pcg(AFun, pct(f), 100*tol, 50);

        % Report a warning if the CG failed to converge
        if (flag ~= 0)
            warning('CHEBFUN:CHEBFUN:inufft2:tolerance',...
                    ['Conjugate gradient method failed to converge to a ' ...
                     'tolerance of %1.2e in %d iterations.  ' ...
                     'Tolerance reached was %1.2e.'],100*tol,50,relres);
        end    
    end

c = plan2(f);
        
p = @plan2;

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
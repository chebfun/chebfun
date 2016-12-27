function [f, p] = nuifft( c, x, omega, type)
% CHEBFUN.NUIFFT   Nonuniform inverse fast Fourier transform
%
% [F, P]= CHEBFUN.NUIFFT( C ) is the same as ifft( C ). C must be a column
% vector. F = P(C) is a planned version of the fast transform. 
%
% F = CHEBFUN.NUIFFT( C, X ) is a nonuniform inverse fast Fourier transform 
% of type 2, which computes A\C, where
% 
%        A_{jk} = exp(-2*pi*1i*X(j)*k/N), 0<=j,K<=N-1. 
%
% C and X must be column vectors of the same length.
%
% F = CHEBFUN.NUIFFT( C, X, 2 ) is the same as CHEBFUN.NUIFFT( C, X )
%
% F = CHEBFUN.NUIFFT( C, OMEGA, 1 ) is a nonuniform inverse fast Fourier 
% transform of type 1, which computes A\C, where
%
%        A_{jk} = exp(-2*pi*1i*j/N*OMEGA(k)), 0<=j,K<=N-1. 
%
% F = CHEBFUN.NUFFT( C, X, 1, TOL), F = CHEBFUN.NUFFT(C, OMEGA, 2, TOL),
% and F = CHEBFUN.NUFFT( C, X, OMEGA, TOL ) are the same as above but with
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
% We do this by solving the normal equations (F'*F)*f=F'*c, where F is the
% NUDFT1 matrix. 

% Plan the NUFFTs:
N = size(omega,1);
[~, p] = chebfun.nufft( c, omega, 1);
[~, p_trans] = chebfun.nufft(c, omega/N, 2); 
p_trans = @(c) conj( p_trans( conj( c ) ) );

% Plan conjugate gradient method on normal equations: 
[f, ~] = pcg(@(c) p_trans( p( c ) ), p_trans(c), 100*tol, 50 );
p = @(c) pcg(@(c) p_ctrans( p_forw(c) ), p_ctrans(c), 100*tol, 50 );
end

function [f,p] = nuifft2( c, x, tol )
% NUIFFT2  Compute the nonuniform IFFT of type 2.
% We do this by solving the normal equations (F'*F)*f=F'*c, where F is the
% NUDFT2 matrix.

% Plan the NUFFTs:
N = size(x,1);
[~, p_forw] = chebfun.nufft(c, x);
[~, p_ctrans] = chebfun.nufft(conj(c), N*x, 1);
p_ctrans = @(y) conj( p_ctrans( conj( y ) ) );

% Plan conjugate gradient method on normal equations: 
[f, ~] = pcg(@(c) p_ctrans( p_forw(c) ), p_ctrans(c), 100*tol, 50 );
p = @(c) pcg(@(c) p_ctrans( p_forw(c) ), p_ctrans(c), 100*tol, 50 );
end
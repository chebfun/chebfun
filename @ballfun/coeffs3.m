function C = coeffs3(f, m, n, p)
% COEFFS3   Trivariate Cheybshev-Fourier-Fourier expansion coefficients of F. 
%   C = COEFFS3(F) returns the tensor of trivariate coefficients.  The
%   coefficients are arranged so they correspond the spherical coordinates
%   radial-azimuthal-polar.
%
%   X = COEFFS3(F, M, N, P) is the same as above but it returns
%   coefficients as an M x N x P tensor.

% Copyright 2018 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( nargin == 1 )
    % Return the coefficients
    C = f.coeffs;
else
    if ( nargin == 2 )
        n = m;
        p = m;
    elseif (nargin == 3 )
        p = n;
    end
    
    %% Copy the diff alias when we solve it below
    
    [mf, nf, pf] = size(f);
    F = f.coeffs;
    
    G = zeros(m,n,pf);

    for k = 1:pf
        G(:,:,k) = chebtech2.alias(trigtech.alias(F(:,:,k).',n).',m);
    end

    G = permute(G,[3,2,1]);
    C = zeros(p,n,m);

    for k = 1:m
       C(:,:,k) = trigtech.alias(G(:,:,k),p); 
    end

    C = permute(C,[3,2,1]);
end
function [coeffs, poles, k] = residue(u, v, k)
%RESIDUE   Partial-fraction expansion (residues).
%   [R, P, K] = RESIDUE(B, A) finds the residues, poles and direct term of a
%   partial fraction expansion of the ratio of two CHEBFUN objects B(s)/A(s).
%   If there are no multiple roots,
%        B(s)       R(1)       R(2)             R(n)
%        ----  =  -------- + -------- + ... + -------- + K(s)
%        A(s)     s - P(1)   s - P(2)         s - P(n)
%   B and A are CHEBFUN objects consisting of a single fun. The residues are
%   returned in the column vector R, the pole locations in column vector P, and
%   the direct terms in the CHEBFUN K. The number of poles is n = length(A) - 1
%   = length(R) = length(P). The direct term CHEBFUN is zero if length(B) <
%   length(A), otherwise length(K) = length(B) - length(A) + 1.
%
%   If P(j) = ... = P(j+m-1) is a pole of multiplicity m, then the expansion
%   includes terms of the form
%                  R(j)        R(j+1)                R(j+m-1)
%                -------- + ------------   + ... + ------------
%                s - P(j)   (s - P(j))^2           (s - P(j))^m
%
%   [B, A] = RESIDUE(R, P, K), with 3 input arguments and 2 output arguments,
%   converts the partial fraction expansion back to the polynomials with CHEBFUN
%   representation in B and A.
%
%   Warning: Numerically, the partial fraction expansion of a ratio of
%   polynomials represents an ill-posed problem. If the denominator polynomial,
%   A(s), is near a polynomial with multiple roots, then small changes in the
%   data, including roundoff errors, can make arbitrarily large changes in the
%   resulting poles and residues. Problem formulations making use of state-space
%   or zero-pole representations are preferable.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( min(size(u)) > 1 || min(size(v)) > 1 )
    error('CHEBFUN:CHEBFUN:residue:quasi', ...
          'Residue does not support CHEBFUN objects with multiple columns.');
end

if ( nargin == 2 )
    
    if ( (numel(u.funs) > 1) || (numel(v.funs) > 1) )
        error('CHEBFUN:CHEBFUN:residue:multipleFuns', ...
              'Residue does not support CHEBFUNs consisting of multiple FUNs.');
    end
    if ( ~isfinite(u) || ~isfinite(v) )
        error('CHEBFUN:CHEBFUN:residue:inf', ...
              'Residue does not support functions with nonzero exponents.');
    end
    
    b = poly(u);
    a = poly(v);
    [coeffs, poles, k] = residue(b, a);
    k = chebfun(@(x) polyval(k, x), max(length(k), 1), 'vectorize');
    
elseif ( nargin == 3 )
    
    % Interpret an empty CHEBFUN as a zero CHEBFUN for the user's convenience.
    if ( isempty(k) )
        k = 0;
    else
        if ( numel(k.funs) > 1 )
            error('CHEBFUN:CHEBFUN:residue:multipleFuns', ...
                  'RESIDUE does not support CHEBFUNs with multiple FUNS.');
        end

        if ( ~isfinite(k) )
            error('CHEBFUN:CHEBFUN:residue:inf', ...
                  'RESIDUE does not support functions which are unbounded.');
        end

        k = poly(k);
    end

    [b, a] = residue(u, v, k);
    coeffs = chebfun(@(x) polyval(b, x), max(length(b), 1), 'vectorize');
    poles = chebfun(@(x) polyval(a, x), max(length(a), 1), 'vectorize');
    
else
    
    error('CHEBFUN:CHEBFUN:residue:args', ...
        'Residue requires either 2 or 3 arguments.');
    
end

end

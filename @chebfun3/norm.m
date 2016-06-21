function [normF, normloc] = norm(f, p)
%NORM   Norm of a CHEBFUN3 object.
%   NORM(F) = NORM(F,'fro') = sqrt(triple integral of abs(F)^2).
%   NORM(F, 1) = NOT IMPLEMENTED AS IT NEEDS SPLITTING CAPABILITIES.
%   NORM(F, inf) = global maximum in modulus.
%   NORM(F, max) is the same as NORM(F, inf).
%   NORM(F, min) = NOT IMPLEMENTED (NEEDS SPLITTING CAPABILITIES).
%
%   The inf norm also returns a second output giving a position where the 
%   max occurs.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( nargin == 1 ) 
    % Default to the Frobenius norm.
    p = 'fro';
end

if ( isempty(f) )  
    % Empty chebfun3.
    normF = [];
    
else
    switch ( p )  % Different cases on different norms.
        case 1
            error('CHEBFUN:CHEBFUN3:norm:norm', ...
                'CHEBFUN3 does not support L1-norm, yet');
        
        case {2}  
            error('CHEBFUN:CHEBFUN3:norm:norm', ...
                'not implemented yet');
    
        case {'fro'}
            % Developer Note: Instead of normF = sqrt(sum3(f.^2)), which 
            % needs explicitly forming the object f.^2, we use HOSVD and 
            % the invariance of the Frobenius norm under unitary modal 
            % multiplication by orthogonal factor quasimatrices. We have: 
            % norm(f, 'fro') = norm(CORE, 'fro') = norm(sv{i}, 'fro') 
            % where CORE is discrete core tensor of the HOSVD of f and
            % sv{i} is the vector of mode-i singular values of f:
            % See chebfun3/hosvd for a reference to the discreteHOSVD paper.
            sv = hosvd(f);
            normF = norm(sv{1}, 'fro');
            
            % Alternatively, a rough approximation could also be computed
            % by sampling f on a discrete tensor F and then computing 
            % Frobenius norm of F.
            
        case {inf, 'inf', 'max'}
            [vals, locs] = minandmax3(f);
            [normF, idx] = max(abs(vals));
            normloc = locs(idx, :);
            
        case {-inf, '-inf', 'min'}
            error('CHEBFUN:CHEBFUN3:norm:norm', ...
                'not implemented.');
            
        case {'op', 'operator'}
            error('CHEBFUN:CHEBFUN3:norm:norm', ...
                'not implemented.');
            
        otherwise
            if ( isnumeric(p) && isreal(p) )
                if ( abs(round(p) - p) < eps )
                    p = round(p); 
                    f = f.^p;
                    if ( ~mod(p,2) )
                        normF = (sum3(f)) .^ (1/p);
                    else
                        error('CHEBFUN:CHEBFUN3:norm:norm', ...
                            'p-norm must have p even for now.');
                    end
                else
                    error('CHEBFUN:CHEBFUN3:norm:norm', ...
                        'CHEBFUN3 does not support this norm.');
                end
            else
                error('CHEBFUN:CHEBFUN3:norm:unknown', 'Unknown norm.');
            end
            
    end
end

end
function [normF, normLoc] = norm(f, n)
%NORM   Norm of a CHEBFUN object.
%   For scalar-valued CHEBFUN objects:
%       NORM(f) = sqrt(integral of abs(f)^2).
%       NORM(f, 2) is the same as NORM(f).
%       NORM(f, 'fro') is also the same as NORM(f).
%       NORM(f, 1) = integral of abs(f).
%       NORM(f, p) = (integral of abs(f)^p)^(1/p).
%       NORM(f, inf) = max(abs(f)).
%       NORM(f, -inf) = min(abs(f)).
%
%   For array-valued CHEBFUN objects:
%       NORM(F) is the Frobenius norm, sqrt(sum(svd(F).^2)).
%       NORM(F, 'fro') is the same as NORM(F).
%       NORM(F, 1) is the maximum of the 1-norms of the columns of F.
%       NORM(F, 2) is the largest singular value of F.
%       NORM(f) = sqrt(integral of abs(f)^2).
%       NORM(F, inf) is the maximum of the 1-norms of the rows of F.
%
% Furthermore, the +\-inf norms for scalar-valued CHEBFUN objects may also
% return a second output, giving the position where the max/min occurs. For
% array-valued CHEBFUN objects, the 1, inf, and p-norms can return as their 2nd
% output the index of the column with the largest norm.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Empty chebfun has norm 0:
if ( isempty(f) )           
    normF = 0;
    return
end

if ( nargin == 1 )
    n = 'fro'; 	% Frobenius norm is the default (2 norm would be much slower)/
end

% Initialise:
normLoc = [];
numCols = size(f.funs{1}, 2);

%% %%%%%%%%%%%%%%%%%%%%%%%%% SCALAR-VALUED CHEBFUNS %%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ( numCols == 1 )
    switch n
        case 1
            if ( nargout == 2 )
                error('CHEBFUN:norm:argout', ...
                        'Cannot return two outputs for 1-norms');
            end
            normF = sum(abs(f));
            
        case {2, 'fro'}
            if ( nargout == 2 )
                error('CHEBFUN:norm:argout', ...
                        'Cannot return two outputs for ''fro''-norms');
            end
            f.isTransposed = 0;
            normF = sqrt(abs(sum(conj(f).*f)));
            
        case {inf, 'inf'}
            if ( isreal(f) )
                [normF, normLoc] = minandmax(f);
                [normF, idx] = max(abs(normF));
                normLoc = normLoc(idx);
            else
                [normF, normLoc] = max(conj(f).*f);
                normF = sqrt(normF);
            end
            
        case {-inf, '-inf'}
            [normF, normLoc] = min(conj(f).*f);
            normF = sqrt(normF);
            
        otherwise
            if ( isnumeric(n) && isreal(n) )
                if ( nargout == 2 )
                    error('CHEBFUN:norm:argout', ...
                            'Cannot return two outputs for p-norms.');
                end
                if ( mod(n, 2) )
                    normF = sum((conj(f).*f).^(n/2))^(1/n);
                else
                    normF = sum(abs(f).^n)^(1/n);
                end
            else
                error('CHEBFUN:norm:unknownNorm', ...
                 'The only matrix norms available are 1, 2, inf, and ''fro''.');
            end
    end
    
%% %%%%%%%%%%%%%%%%%%%%%%%%%% ARRAY-VALUED CHEBFUNS %%%%%%%%%%%%%%%%%%%%%%%%%%%%
else 
    
    switch n
        case 1
            f = mat2cell(f);
            normF = zeros(1, numel(f));
            for k = 1:numel(f)
                normF(k) = norm(f{k}, 1);
            end
            [normF, normLoc] = max(normF);
            
        case 2
            if (nargout == 2 )
                error('CHEBFUN:norm:argout', ...
                    'Cannot return two outputs for quasimatrix 2-norms.');
            end
            s = svd(f, 0);
            normF = s(1);
            
        case 'fro'
            if ( nargout == 2 )
                error('CHEBFUN:norm:argout', ...
                    'Cannot return two outputs for quasimatrix ''fro''-norms.');
            end
            % Find integration dimension: 1 if column, 2 if row:
            f.isTransposed = 0;
            normF = sqrt(sum(sum(f.*conj(f))));
            
        case {'inf', inf}
            f.isTransposed = 0;
            [normF, normLoc] = max(sum(abs(f), 2));
            
        case {'-inf', -inf}
            f.isTransposed = 0;
            [normF, normLoc] = min(sum(abs(f), 2));
            
        otherwise
            if ( isnumeric(n) && isreal(n) )
                f.isTransposed = 0;
                [normF, normLoc] = max(sum(abs(f).^n, 2));
                normF = normF^(1/n);
            else
                error('CHEBFUN:norm:unknownNorm', ...
                 'The only matrix norms available are 1, 2, inf, and ''fro''.');
            end
    end
    
end

% Discard possible imaginary rounding errors:
normF = real(normF);

end

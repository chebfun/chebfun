function [normF, normLoc] = norm(F, n)
%NORM   Norm of a chebfun object.
%   For scalar-valued chebfun objects:
%       NORM(f) = sqrt(integral of abs(f)^2).
%       NORM(f, 2) is the same as NORM(f).
%       NORM(f, 'fro') is also the same as NORM(f).
%       NORM(f, 1) = integral of abs(f).
%       NORM(f, inf) = max(abs(f)).
%       NORM(f, -inf) = min(abs(f)).
%
%   For array-valued chebfun objects:
%       NORM(F) is the Frobenius norm, sqrt(sum(svd(F).^2)).
%       NORM(F, 1) is the maximum of the 1-norms of the columns of F.
%       NORM(F, 2) is the largest singular value of F.
%       NORM(F, inf) is the maximum of the 1-norms of the rows of F.
%       NORM(F, 'fro') is the same as NORM(F).
%
% Furthermore, the +\-inf norms for scalar-vaued chebfun objects may also return
% a second input, giving the position where the max/min occurs. For array-valued
% chebfun objects, the 1, inf, and p-norms can return as their 2nd output the
% index of the column with the largest norm.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( nargin == 1 )
    n = 'fro'; 	% Frobenius norm is the default
end             % (2 norm would be much slower)


normLoc = [];

if ( isempty(F) )           % Empty chebfun has norm 0.
    normF = 0;
elseif ( size(F, 2) == 1 )  % A is a scalar-valued chebfun.
    switch n
        case 1
            if ( nargout == 2 )
                error('CHEBFUN:norm:argout',...
                        'Cannot return two outputs for 1-norms');
            end
            absA = abs(F);
            normF = sum(absA);
        case {2, 'fro'}
            if ( nargout == 2 )
                error('CHEBFUN:norm:argout',...
                        'Cannot return two outputs for ''fro''-norms');
            end
            if ( F.isTransposed )
                normF = sqrt(abs(F*F'));
            else
                normF = sqrt(abs(F'*F));
            end
        case {inf, 'inf'}
            if ( isreal(F) )
                [normF, normLoc] = minandmax(F);
                [normF, idx] = max([-normF(1), normF(2)]);
                normLoc = normLoc(idx);
            else
                [normF, normLoc] = max(conj(F).*F);
                normF = sqrt(normF);
            end
        case {-inf, '-inf'}
            [normF, normLoc] = min(abs(F));
        otherwise
            if ( isnumeric(n) && isreal(n) )
                if ( nargout == 2 )
                    error('CHEBFUN:norm:argout',...
                            'Cannot return two outputs for p-norms');
                end
                if ( mod(n, 2) )
                    normF = sum((conj(F).*F).^(n/2))^(1/n);
                else
                    normF = sum(abs(F).^n)^(1/n);
                end
            else
                error('CHEBFUN:norm:unknown', 'Unknown norm.');
            end
    end

else                        % A is a vector-valued:
    % [TODO]: Implement norms of vector-valued chebfun objects.
    error('Norm not yet implemented')
    switch n
        case 1
            normF = zeros(numel(F));
            for k = 1:numel(F)
                normF(k) = norm(F(k), 1);
            end
            [normF, normLoc] = max(norm(F));
        case 2
            if nargout == 2
                error('CHEBFUN:norm:argout',...
                    'Cannot return two outputs for quasimatrix 2-norms');
            end
            s = svd(F,0);
            normF = s(1);
        case 'fro'
            % Find integration dimension: 1 if column, 2 if row
            if nargout == 2
                error('CHEBFUN:norm:argout', ...
                    'Cannot return two outputs for quasimatrix ''fro''-norms');
            end
            dim = 1 + double(F(1).trans);
            normF = sqrt( sum( sum(F.*conj(F),dim) ) );
        case {'inf',inf}
            [normF, normLoc] = max(sum(abs(F),2));
        case {'-inf',-inf}
            [normF, normLoc] = min(sum(abs(F),2));
        otherwise
            if isnumeric(n) && isreal(n)
%                 normA = max(sum((conj(A).*A).^(n/2)))^(1/n);
%                 normA = max(sum(abs(A).^n))^(1/n);
                [normF, normLoc] = max(sum(abs(F).^n));
                normF = normF^(1/n);
            else
                error('CHEBFUN:norm:unknown2','Unknown norm');
            end
    end
end

% Discard possible imaginary rounding errors:
normF = real(normF);

end

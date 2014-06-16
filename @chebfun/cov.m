function out = cov(f, g, varargin)
%COV   Covariance of a CHEBFUN.
%   COV(F) is the same as VAR(F) if F is a scalar-valued CHEBFUN.
%   COV(F) returns the covariance of the array-valued CHEBFUN F. 
%   COV(F, G) returns the covariance matrix of the columns of F and G.
%
% See also VAR, MEAN.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Empty case:
if ( isempty(f) )
    out = NaN;
    return
end

% Error checking:
if ( (nargin == 2) && ~all(size(f) == size(g)) )
    error('CHEBFUN:CHEBFUN:cov:size',' CHEBFUN dimensions do not agree.');
end
if ( nargin == 3 )
    error('CHEBFUN:CHEBFUN:cov:nargin', ...
        'CHEBFUN/COV does not support normalization.');
end

% Deal with row CHEBFUN objects:
if ( f(1).isTransposed )
    if ( nargin == 1 )
        out = transpose(cov(transpose(f)));
    else
        out = transpose(cov(transpose(f), transpose(g)));
    end
    return
end

% Conditional on COV(f) and COV(f, g).
if ( nargin == 1 ) % COV(f)
    
    if ( numColumns(f) == 1 )
        % The covariance of a scalar-valued CHEBFUN is the same as the variance:
        out = var(f);
        return
        
    else
        % Array-valued CHEBFUN or quasimatrix.
        
        Y = f - mean(f);
        out = diag(mean(Y.*conj(Y)));
        % Convert Y to a cell array of scalar-valued CHEBFUN objects.
        Y = mat2cell(Y);
        % Loop over each of the columns:
        for j = 1:numel(Y)
            for k = j+1:numel(Y)
                % Compute the scaled inner product of the jth and kth columns:
                out(j,k) = mean(Y{j}.*conj(Y{k}));
                % Use symmetry:
                out(k,j) = conj(out(j,k));
            end
        end
        
    end
        
else               % COV(f, g)
    
    % Convert to cell arrays of scalar-valued CHEBFUN objects:
    Y = cheb2cell(f - mean(f));
    Z = cheb2cell(g - mean(g));
    % Initialise output matrix:
    out = zeros(numel(f));
    % Loop over each of the columns:
    for j = 1:numel(Y)
        for k = 1:numel(Y)
            out(j,k) = mean(Y{j}.*conj(Z{k}));
        end
    end

end

end

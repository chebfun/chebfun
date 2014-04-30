function [normA, normLoc] = norm(A, n)
%NORM   Norm of a CHEBMATRIX object.
%   NORM(A) computes the norm of the CHEBMATRIX object A.
%
%   If A contains only CHEBFUN and DOUBLE objects, A is converted to a
%   QUASIMATRIX, and CHEBFUN/NORM is called.
%
%   If not [TODO].
%
%   See also CHEBMATRIX, CHEBFUN/NORM.

%  Copyright 2014 by The University of Oxford and The Chebfun Developers.
%  See http://www.chebfun.org for Chebfun information.

% Empty CHEBMATRIX has norm 0:
if ( isempty(A) )
    normA = 0;
    return
end

if ( nargin == 1 )
    n = 'fro'; 	% Frobenius norm is the default.
end

% Initialise:
normLoc = [];
temp = 1;
sz = size(A, 1)*size(A, 2);

% Check if A contains only CHEBFUN and DOUBLE objects.
for j = 1:sz
    if  ( isa(A.blocks{j}, 'chebfun') | isa(A.blocks{j}, 'double') )
    else
        temp = 0;    
    end
end

% If so, convert A to a QUASIMATRIX.
if temp == 1
    A = chebfun(A);
    numCols = size(A, 2);
    % If A is actually simply a CHEBFUN, that is a QUASIMATRIX with 
    % one column.
    if ( numCols == 1 )
        switch n
            case {1, 2, 'fro'}
            normA = norm(A, n);
        
            case {inf, 'inf', -inf, '-inf'} 
            [normA, normLoc] = norm(A, n); 
            
            otherwise
                if ( isnumeric(n) && isreal(n) )
                    [normA, normLoc] = norm(A, n); 
                else
                error('CHEBMATRIX:norm:unknownNorm', ...
                 'The only matrix norms available are 1, 2, inf, and ''fro''.');
            end
        end
    % If A is a QUASIMATRIX with more that one column.
    else
        switch n
            case {2, 'fro'}
            normA = norm(A, n);
        
            case {1, inf, 'inf', -inf, '-inf'} 
            [normA, normLoc] = norm(A, n);
            
            otherwise
                if ( isnumeric(n) && isreal(n) )
                    [normA, normLoc] = norm(A, n); 
                else
                error('CHEBMATRIX:norm:unknownNorm', ...
                 'The only matrix norms available are 1, 2, inf, and ''fro''.');
            end
        end
    end
   
% If not [TODO].
else
    normA = 0;
    normLoc = 0;
end

end
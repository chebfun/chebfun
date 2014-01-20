function I = identity(A)
%IDENTITY  Identity chebmatrix following a given variable structure.
%   Suppose A is a square chebmatrix (identical row and column block sizes,
%   so that A*A is well defined). Then I = IDENTITY(A) returns an identity
%   operator such that I*A = A*I = A for all identically sized A. The
%   diagonal blocks of A are all either identity blocks or scalar 1's. 
%
%   See also CHEBMATRIX, CHEBMATRIX.BLOCKSIZES.
    
%  Copyright 2013 by The University of Oxford and The Chebfun Developers.
%  See http://www.chebfun.org for Chebfun information.

% This method is not called chebmatrix.eye because then the static method
% linop.eye would overload the desired chebmatrix behavior. 

if ( size(A,1) ~= size(A,2) )
    error('Chebmatrix must be square.')
end

I = chebmatrix(A.blocks);
isFun = isFunVariable(A);
n = size(A,1);
d = A.domain;
for i = 1:n
    
    % Diagonal block maps variable type to itself--just two cases.
    if ( isFun(i) )
        I.blocks{i,i} = linBlock.eye(d);
    else
        I.blocks{i,i} = 1;
    end
    
    % Off diagonal blocks have 2x2 cases.
    for j = [1:i-1 i+1:n]
        if ( isFun(i) )
            if ( isFun(j) )
                % fun -> fun
                I.blocks{i,j} = linBlock.zeros(d);
            else
                % scalar -> fun
                I.blocks{i,j} = chebfun(0,d);
            end
        else
            if ( isFun(j) )
                % fun -> scalar
                I.blocks{i,j} = linBlock.zero(d);
            else
                % scalar -> scalar
                I.blocks{i,j} = 0;
            end
        end
    end
    
end
            
end

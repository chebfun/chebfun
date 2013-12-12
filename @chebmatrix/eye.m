function I = blockEye(A)
%EYE    Identity chebmatrix following a given structure.

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

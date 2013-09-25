function A = matrix(L,dim,varargin)
% MATRIX(A,DIM) returns a discretization of the chebmatrix A based on dimension DIM.
%
% MATRIX(A,DIM,TYPE) uses the discretization type TYPE.

A = matrixBlocks(L,dim,varargin{:});

try
    A = cell2mat(A);
catch
    %TODO: If this is checked using the cat methods, should this even be here?
    error('Block sizes are not compatible.')
end

end

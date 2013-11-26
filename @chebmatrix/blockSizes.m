function varargout = blockSizes(A)
%BLOCKSIZES Sizes of the blocks within the chebmatrix.
% BLOCKSIZES(L) returns a cell of 1x2 size vectors.
% [M, N] = BLOCKSIZES(A) returns two matrices of row/column sizes.

if ( nargout <= 1 )
    varargout = {cellfun(@size, A.blocks, 'uniform', false)};
else
    varargout{1} = cellfun(@(x) size(x, 1), A.blocks);
    varargout{2} = cellfun(@(x) size(x, 2), A.blocks);
end

end

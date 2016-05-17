function varargout = blockSizes(A)
%BLOCKSIZES Sizes of the blocks within a chebmatrix.
%   BLOCKSIZES(L) returns a cell of 1x2 size vectors. Each entry is one of
%   these:
%     [Inf,Inf] : operator block (maps function to function)
%     [  1,Inf] : functional block (maps function to scalar)
%     [Inf,  1] : chebfun block (maps scalar to function)
%     [  1,  1] : scalar block (maps scalar to scalar)
%
%   [M, N] = BLOCKSIZES(A) returns two matrices of row/column sizes.
%
% See also CHEBMATRIX, CHEBMATRIX.SIZE.
    
% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( nargout <= 1 )
    varargout = {cellfun(@size, A.blocks, 'uniform', false)};
else
    varargout{1} = cellfun(@(x) size(x, 1), A.blocks);
    varargout{2} = cellfun(@(x) size(x, 2), A.blocks);
end

end

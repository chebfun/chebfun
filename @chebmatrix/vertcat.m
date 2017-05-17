function C = vertcat(varargin)
%VERTCAT   Horizontally concatenate chebmatrices.
%   Z = [A; B; C; ...] vertically combines the chebmatrices, operator blocks,
%   chebfuns, and scalars given in the call, if their column sizes are
%   compatible.
%
% See also CHEBMATRIX.CAT, CHEBMATRIX.HORZCAT.

%  Copyright 2017 by The University of Oxford and The Chebfun Developers.
%  See http://www.chebfun.org/ for Chebfun information.

% Remove empty arguments.
isemp = cellfun(@isempty,varargin);
varargin(isemp) = [];
nargin = length(varargin);

% First, process the blocks of the inputs, disregarding sizes. This
% simplifies the logic. The entries of the results are cells holding the
% relevant blocks to be concatenated.
blocks = cell(1, nargin);
for j = 1:nargin
    item = varargin{j};
    if ( isa(item, 'chebmatrix') )
        blocks{j} = item.blocks;
    elseif ( isa(item, 'linBlock') )
        blocks{j} = { item };
    elseif ( ~isempty(item) )
        blocks{j} = num2cell(item);
    end
end

% Now check the column sizes. Create an empty vessel if OK.
[m, n] = cellfun(@size, blocks);
if ( all(n==n(1)) )
    B = cell( sum(m), n(1) );
else
    error('CHEBFUN:CHEBMATRIX:vertcat:sizeMismatch', ...
        'Incompatible column sizes.')
end

% Now we need to flatten out all the inner nested cell divisions, leaving
% just a cell of blocks. 
cs = [0 cumsum(m)];
for j = 1:nargin
    B(cs(j)+(1:m(j)),:) = blocks{j};
end

% This step will perform domain compatibility checking.
C = chebmatrix( B );

% Finally, check block size compatibility.
[row, col] = blockSizes(C);
if ( any( any( bsxfun(@ne, col(1,:), col) ) ) )
         error('CHEBFUN:CHEBMATRIX:vertcat:blockSizeMismatch', ...
            'Block column sizes must be the same.')
end

end

function C = horzcat(varargin)
%HORZCAT   Horizontally concatenate chebmatrices.
%   Z = [A,B,C,...] horizontally combines the chebmatrices, operator
%   blocks, chebfuns, and scalars given in the call, if their row sizes
%   are compatible. 
%
% See also CHEBMATRIX.CAT, CHEBMATRIX.VERTCAT.

%  Copyright 2016 by The University of Oxford and The Chebfun Developers.
%  See http://www.chebfun.org/ for Chebfun information.

% Remove empty arguments.
isemp = cellfun(@isempty,varargin);
varargin(isemp) = [];
nargin = length(varargin);

% First, process the blocks of the inputs, disregarding sizes. This
% simplifies the logic. The entries of the results are cells holding the
% relevant blocks to be concatenated.
blocks = cell(1, nargin);
for n = 1:nargin
    item = varargin{n};
    if ( isa(item, 'chebmatrix') )   % already a cell
        blocks{n} = item.blocks;
    elseif ( isa(item, 'linBlock') || isa(item,'chebfun') )
        blocks{n} = { item };        % encase in a cell
    elseif ( ~isempty(item) )   
        blocks{n} = num2cell(item);  % encase in cell(s)
    end
end

% Now check the row sizes. Create an empty vessel if OK.
[m,n] = cellfun(@size,blocks);
if ( all(m==m(1)) )
    B = cell( m(1), sum(n) );
else
    error('CHEBFUN:CHEBMATRIX:horzcat:sizeMismatch', ...
        'Incompatible row sizes.')
end

% Now we need to flatten out all the inner nested cell divisions, leaving
% just a cell of blocks. 
cs = [0 cumsum(n)];
for j = 1:nargin
    B(:,cs(j)+(1:n(j))) = blocks{j};
end

% This step will perform domain compatibility checking.
C = chebmatrix( B );

% Finally, check block size compatibility.
[row, col] = blockSizes(C);
if ( any( any( bsxfun(@ne, row(:,1), row) ) ) )
    error('CHEBFUN:CHEBMATRIX:horzcat:blockSizeMismatch', ...
        'Block row sizes must be the same.')
end

end


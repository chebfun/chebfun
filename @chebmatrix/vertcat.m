function C = vertcat(varargin)

%  Copyright 2013 by The University of Oxford and The Chebfun Developers.
%  See http://www.chebfun.org for Chebfun information.

% Remove empty arguments.
isemp = cellfun(@isempty,varargin);
varargin(isemp) = [];
nargin = length(varargin);

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

[m,n] = cellfun(@size,blocks);

if ( all(n==n(1)) )
    B = cell( sum(m), n(1) );
else
    error('Incompatible column sizes.')
end

cs = [0 cumsum(m)];
for j = 1:nargin
    B(cs(j)+(1:m(j)),:) = blocks{j};
end

% This step will include checking domain compatibility.
C = chebmatrix( B );

% Check block size compatibility.
[row, col] = blockSizes(C);
if ( any( any( bsxfun(@ne, col(1,:), col) ) ) )
         error('Block column sizes must be the same.')
end

end

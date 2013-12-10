function C = horzcat(varargin)

%  Copyright 2013 by The University of Oxford and The Chebfun Developers.
%  See http://www.chebfun.org for Chebfun information.

% Remove empty arguments.
isemp = cellfun(@isempty,varargin);
varargin(isemp) = [];
nargin = length(varargin);

blocks = cell(1, nargin);
for n = 1:nargin
    item = varargin{n};
    if ( isa(item, 'chebmatrix') )
        blocks{n} = item.blocks;
    elseif ( isa(item, 'linBlock') )
        blocks{n} = { item };
    elseif ( ~isempty(item) )
        blocks{n} = num2cell(item);
    end
end

[m,n] = cellfun(@size,blocks);

if ( all(m==m(1)) )
    B = cell( m(1), sum(n) );
else
    error('Incompatible row sizes.')
end

cs = [0 cumsum(n)];
for j = 1:nargin
    B(:,cs(j)+(1:n(j))) = blocks{j};
end

% This step will include checking domain compatibility.
C = chebmatrix( B );

% Check block size compatibility.
[row, col] = blockSizes(C);
if ( any( any( bsxfun(@ne, row(:, 1), row) ) ) )
         error('Block row sizes must be the same.')
end

end


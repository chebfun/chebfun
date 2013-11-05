function C = cat(dim, varargin)

% Any singleton operator blocks must encased in cells to match up with
% chebmatrix blocks.
blocks = cell(1, nargin-1);
for n = 1:nargin-1
    if ( isa(varargin{n}, 'chebmatrix') )
        blocks{n} = varargin{n}.blocks;
    elseif ( isa(varargin{n}, 'linBlock') )
        blocks{n} = varargin(n);
    elseif ( ~isempty(varargin{n}) )
        blocks{n} = num2cell(varargin{n});
    end
end

% This step will include checking domain compatibility.
C = chebmatrix( cat(dim, blocks{:}) );

% Check block size compatibility.
[row, col] = blockSizes(C);

if ( dim == 1 )
    if ( any( any( bsxfun(@ne, col(1, :), col) ) ) )
        error('Block column sizes must be the same.')
    end
elseif ( dim == 2 )
    if ( any( any( bsxfun(@ne, row(:, 1), row) ) ) )
        error('Block row sizes must be the same.')
    end
else
    error('Must concatenate along dimension 1 or 2.')
end  

end

function A = matrix(L,varargin)

dsc = L.discretizationType(L,varargin{:});

A = cell(size(L.blocks));
for i = 1:numel(A)
    item = L.blocks{i};
    if isa(item,'linBlock')
        dsc.source = item;
        A{i} = matrix(dsc);
    elseif isa(item,'chebfun')
        A{i} = dsc.toValues(item);
        if ( item.isTransposed )
            A{i} = A{i}.';
        end
    elseif isnumeric(item)
        A{i} = item;
    else
        error('Unrecognized block type.')
    end
end

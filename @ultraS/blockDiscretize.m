function [L, S] = blockDiscretize(disc, block)
%  Copyright 2013 by The University of Oxford and The Chebfun Developers.
%  See http://www.chebfun.org for Chebfun information.
S = [];
if (isa(block, 'operatorBlock') )
    if ( ~isempty(disc.coeffs) )
        [L, S] = quasi2USdiffmat(disc);
    else
        error
    end
elseif ( isa(block, 'functionalBlock') )
    dim = disc.dimension;
    dom = disc.domain;
    collocDisc = colloc2(block,dim,dom);
    L = matrix(collocDisc);
    cumsumDim = [0, cumsum(dim)];
    numInts = numel(dom) - 1;
    tmp = cell(1, numInts);
    for k = 1:numInts
        Lk = L(cumsumDim(k) + (1:dim(k)));
        tmp{k} = flipud(chebtech2.coeffs2vals(Lk.')).';
    end
    L = cell2mat(tmp);
elseif ( isa(block, 'chebfun') )
    L = toValues(disc, block);
    if ( block.isTransposed )
        error % TODO: ?
        L = L.';
    end
elseif ( isnumeric(block) )
    L = block;
else
    error('Unrecognized block type.')
end
end

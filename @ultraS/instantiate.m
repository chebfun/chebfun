function [L, S] = instantiate(disc,data)
%INSTANTIATE Convert an item to discrete form in ULTRAS.
%   [L,S] = INSTANTIATE(DISC,DATA) converts each item DATA{k} to discrete form using the
%   information in discretization DISC. The result M is a cell array.
%
%   Each item may be:
%      linBlock (becomes a matrix)
%      chebfun (becomes a vector)
%      numeric (not changed)

%  Copyright 2013 by The University of Oxford and The Chebfun Developers.
%  See http://www.chebfun.org for Chebfun information.

if ( iscell(data) )
    [L,S] = cellfun(@(x) instantiateOne(disc,x),data,'uniform',false);
else
    [L,S] = instantiateOne(disc,data);
end

end

function [L,S] = instantiateOne(disc,item)

if (isa(item, 'operatorBlock') )
    if ( ~isempty(disc.coeffs) )
        [L, S] = quasi2USdiffmat(disc);
    else
        error('CHEBFUN:ultraS:fail', ...
            'ultraS cannot represent this operator. Suggest you use colloc2.')   
    end 
elseif ( isa(item, 'functionalBlock') )
    dim = disc.dimension;
    dom = disc.domain;
    collocDisc = colloc2(item,dim,dom);
    L = matrix(collocDisc);
    cumsumDim = [0, cumsum(dim)];
    numInts = numel(dom) - 1;
    tmp = cell(1, numInts);
    for k = 1:numInts
        Lk = L(cumsumDim(k) + (1:dim(k)));
        tmp{k} = flipud(chebtech2.coeffs2vals(Lk.')).';
    end
    L = cell2mat(tmp);
    S = zeros(size(L));
elseif ( isa(item, 'chebfun') )
    L = toValues(disc, item);
    if ( item.isTransposed )
        error % TODO: ?
        L = L.';
    end
    S = zeros(size(L));
elseif ( isnumeric(item) )
    L = item;
    S = 1;
else
    error('Unrecognized item type.')
end
end

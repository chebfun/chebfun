function [L, S] = instantiate(disc, data)
%INSTANTIATE Convert an item to discrete form in ULTRAS.
%   [L, S] = INSTANTIATE(DISC, DATA) converts each item DATA{k} to discrete form
%   using the information in discretization DISC. The result M is a cell array.
%
%   Each item may be:
%      linBlock (becomes a matrix)
%      chebfun (becomes a vector)
%      numeric (not changed)

%  Copyright 2014 by The University of Oxford and The Chebfun Developers.
%  See http://www.chebfun.org for Chebfun information.

% TODO: Why two outputs L and S? What do they mean? You say M is a cell-array,
% this seems to be oudated info. Aren't L and S just matrices?

if ( iscell(data) )
    % Loop through cells
    [L, S] = cellfun(@(x) instantiateOne(disc, x), data, 'uniform', false);
else
    [L, S] = instantiateOne(disc,data);
end

end

function [L,S] = instantiateOne(disc, item)
% Instantiate one block of data.
%
% TODO: Why two outputs L and S? What do they mean?

if (isa(item, 'operatorBlock') )
    % Convert a square block
    if ( ~isempty(disc.coeffs) )
        % Coefficients of the block are available, convert to a diffmat.
        [L, S] = quasi2USdiffmat(disc);
    else
        error('CHEBFUN:ultraS:fail', ...
            'ultraS cannot represent this operator. Suggest you use colloc2.')   
    end 
elseif ( isa(item, 'functionalBlock') )
    % Convert a row block
    
    % TODO: More documentation please. AB, 12/2/14
    dim = disc.dimension;
    dom = disc.domain;
    collocDisc = colloc2(item, dim, dom);
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
    % Block is a CHEBFUN. Convert to value space.
    L = toValues(disc, item);
    if ( item.isTransposed )
        error % TODO: ?
        L = L.';
    end
    S = zeros(size(L));
elseif ( isnumeric(item) )
    % Block is numeric, don't need to do much.
    L = item;
    S = 1;
else
    error('Unrecognized item type.')
end
end

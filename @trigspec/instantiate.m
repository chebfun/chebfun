function [M, S] = instantiate(disc)
%INSTANTIATE   Convert a TRIGSPEC discretization to discrete form.
%   M = INSTANTIATE(DISC) converts each item DISC.SOURCE to discrete form
%   using the information in discretization DISC. The result M is return a cell
%   array if DISC.SOURCE has more than one component.
%
%   [M, S] = INSTANTIATE(DISC) retusn a second output, S, which is an empty
%   cell-array.
%
%   DISC.SOURCE may be one or a cell array of:
%      linBlock (becomes a matrix)
%      chebfun (becomes a vector)
%      numeric (not changed)

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

data = disc.source;
if ( isa(data, 'chebmatrix') )
    data = data.blocks;
end

if ( iscell(data) )
    M = cell(size(data));
    S = cell(size(data));
    for j = 1:size(data, 1)
        for k = 1:size(data, 2)
            discJK = extractBlock(disc, j, k);
            [M{j,k}, S{j,k}] = instantiate(discJK);
        end
    end
    return
else
    [M, S] = instantiateOne(disc, data);
end

end

function [M, S] = instantiateOne(disc, item)
% Instantiate one block of data.

if ( isa(item, 'operatorBlock') )
    % Convert a square block.
    
    if ( ~isempty(disc.coeffs) )
        % Coefficients of the block are available, convert to a diffmat.
        [M, S] = quasi2diffmat(disc);
    else
        error('CHEBFUN:TRIGSPEC:instantiate:fail', ...
         'TRIGSPEC cannot represent this operator. Suggest you use TRIGCOLLOC.')
    end
    
elseif ( isa(item, 'chebfun') )
    % Block is a CHEBFUN. Convert to value space.
    
    M = toValues(disc, item);
    if ( item.isTransposed )
        M = M.';
    end
    S = [];
    
elseif ( isnumeric(item) )
    % Block is numeric, don't need to do much.
    
    M = item;
    S = [];
    
else
    
    error('CHEBFUN:TRIGSPEC:instantiate:inputType', ...
        'Unrecognized item type.')
    
end

end

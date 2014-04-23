function M = instantiate(disc, data)
%INSTANTIATE Convert an item to discrete form in COLLOC2.
%   M = INSTANTIATE(DISC, DATA) converts each item DATA{k} to discrete form
%   using the information in discretization DISC. The result M is a cell array.
%
%   Each item may be:
%      linBlock (becomes a matrix)
%      chebfun (becomes a vector)
%      numeric (not changed)
%
%   See also: COLLOC2/MATRIX

%  Copyright 2013 by The University of Oxford and The Chebfun Developers.
%  See http://www.chebfun.org for Chebfun information.

if ( nargin < 2 )
    data = disc.source.blocks;
end

if ( iscell(data) )
    M = cell(size(data));    
    for j = 1:size(data, 1)
        for k = 1:size(data, 2)
            discJK = extractBlock(disc, j, k);
            M{j,k} = instantiate(discJK);
        end
    end
else
    M = instantiateOne(data);
end


    function A = instantiateOne(item)
    % Instantiate individual component of a discretization+chebmatrix combo.
        if ( isa(item, 'linBlock') )
            disc.source = item;
            A = disc.source.stack( disc );
        elseif ( isa(item, 'chebfun') )
            A = disc.toValues(item);
        elseif ( isnumeric(item) )
            A = item;
        else
            error('CHEBFUN:COLLOC2:instantiate:instantiateOne', 'Unrecognized item.')
        end
    end

end
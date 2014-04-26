function [M, S] = instantiate(disc)
%INSTANTIATE   Convert a COLLOC discretization to discrete form.
%   M = INSTANTIATE(DISC) converts each item DISC.SOURCE to discrete form
%   using the information in discretization DISC. The result M is a cell
%   array if DISC.SOURCE has more than one component.
%
%   DISC.SOURCE may be one or a cell array of:
%      linBlock (becomes a matrix)
%      chebfun (becomes a vector)
%      numeric (not changed)
%
% See also: MATRIX

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org for Chebfun information.

% TODO: Document S.

data = disc.source;
if ( isa(data, 'chebmatrix') )
    data = data.blocks;
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
    M = instantiateOne(disc, data);
end

S = {};

end

function A = instantiateOne(disc, item)
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

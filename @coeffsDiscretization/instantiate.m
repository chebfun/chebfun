function [M, S] = instantiate(disc)
%INSTANTIATE   Convert a COEFFSDISCRETIZATION to discrete form.
%   M = INSTANTIATE(DISC) converts each item DISC.SOURCE to discrete form
%   using the information in discretization DISC. The result M is return a cell
%   array if DISC.SOURCE has more than one component.
%
%   [M, S] = INSTANTIATE(DISC) retusn a second output, S, which is a cell array
%   containing the dscrete form of the conversion operator for each block
%   of DISC.SOURCE. Only used by ULTRAS class.
%
%   DISC.SOURCE may be one or a cell array of:
%      linBlock (becomes a matrix)
%      chebfun (becomes a vector)
%      numeric (not changed)

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
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
        [M, S] = convertOperator(disc, item);
    end
    
elseif ( isa(item, 'functionalBlock') )
    % Convert a row block.
    
    [M, S] = convertFunctional(disc, item);
    
elseif ( isa(item, 'chebfun') )
    % Block is a CHEBFUN. Convert to value space.
    
    M = toValues(disc, item);
    if ( item.isTransposed )
        M = M.';
    end
    S = zeros(size(M));
    
elseif ( isnumeric(item) )
    % Block is numeric, don't need to do much.
    
    M = item;
    S = 1;
    
else
    
    error('COEFFSDISCRETIZATION:instantiate:inputType', ...
        'Unrecognized item type.')
    
end

end

function [M, S] = convertFunctional(disc, item)
    % Developer note: In general we can't represent functional
    % blocks via coeffs. To get around this we instantiate a
    % VALSDISCRETIZAION and convert it to coefficient space
    % using COEFFS2VALS(). (Note it's COEFFS2VALS() rather than
    % VALS2COEFFS() because it's a right-multiply (I think..).)

    % For convenience:
    dim = disc.dimension;
    dom = disc.domain;

    % Create a TECH and a VALSDISCRETIZATION:
    tech = disc.returnTech;
    tech = tech();
    valsDisc = tech.returnValsDisc;
    valsDisc = valsDisc(item, dim, dom);
    M = matrix(valsDisc);

    % Convert from value-space to coeff-space using COEFFS2VALS.
    cumsumDim = [0, cumsum(dim)];
    tmp = cell(1, numel(dom)-1);
    for l = 1:numel(tmp)
        Ml = M(cumsumDim(l) + (1:dim(l)));
        Ml = rot90(Ml);
        tmp{l} = rot90(tech.coeffs2vals(Ml), -1);
    end

    M = cell2mat(tmp);
    S = zeros(size(M));
end
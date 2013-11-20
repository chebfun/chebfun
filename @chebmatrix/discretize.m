function A = discretize(L, disc)

% Construct a single matrix based on DIM-dimensional blocks.

if ( isnumeric(disc) )
    dim = disc;
    % TODO: Choose a proper default disc.
    disc = colloc2(L);
    disc.dimension = dim;
end

disc.domain = L.domain;
if ( numel(disc.dimension) == 1 )
    disc.dimension = repmat(disc.dimension, 1, numel(disc.domain)-1);
end

A = discretize(disc);

end

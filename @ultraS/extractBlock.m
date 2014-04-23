function disc = extractBlock(disc, j, k)
    disc.source.blocks = disc.source.blocks{j,k};
    disc.coeffs = disc.coeffs{j, k};
    disc.outputSpace = disc.outputSpace(j);
    if ( numel(disc.inputDimension) > 1 )
        disc.inputDimension = disc.inputDimension(j,k);
    end
    disc.dimension = disc.dimension + disc.inputDimension;
end

function disc = extractBlock(disc, j, k)
    disc.source.blocks = disc.source.blocks{j,k};
    disc.inputDimension = disc.inputDimension(j,k);
    disc.dimension = disc.dimension + disc.inputDimension;
end
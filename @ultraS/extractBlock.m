function disc = extractBlock(disc, j, k)
    disc.source.blocks = disc.source.blocks{j,k};
    disc.coeffs = disc.coeffs{j, k};
    disc.outputSpace = disc.outputSpace(j);
    disc.inputDimension = disc.inputDimension(j,k);
    disc.dimension = disc.dimension + disc.inputDimension;
end

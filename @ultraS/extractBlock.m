function disc = extractBlock(disc, j, k)
    disc.source.blocks = disc.source.blocks{j,k};
    disc.coeffs = disc.coeffs{j, k};
    disc.outputSpace = disc.outputSpace(j);
    if ( numel(disc.inputDimensionAdjustment) > 1 )
        disc.inputDimensionAdjustment = disc.inputDimensionAdjustment(j,k);
    end
    disc.dimension = disc.dimension + disc.inputDimensionAdjustment;
end

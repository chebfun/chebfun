function ok = validate(disc)

ok = true;
if isempty(disc.dimension)
    error('Dimension of discretization not specified.')
end
if ( length(disc.dimension) ~= disc.numIntervals )
    error('Discretization lengths improperly specified for subintervals.')
end

end

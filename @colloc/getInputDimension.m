function space = getInputDimension(L)

if ( isa(L, 'linop') )
    space = max(getDiffOrder(L), [], 1);
    space = max(space, 0);
    space = repmat(space, size(L, 1), 1);
else
    space = zeros(size(L));
end

end
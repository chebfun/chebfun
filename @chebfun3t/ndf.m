function out = ndf(F)
%NDF   Number of degrees of freedom (parameters) needed to represent a Chebfun3t.

if ( isempty(F) )
    out = 0;
else
    out = prod(size(F.coeffs));
end

end

% F.rank = number of terms in the slice decomposition.
% length(F.rows(:,i)) = the number of entries in each vector that
% corresponds to a row of F.
% rank(F.slices{i}) = rank of the i-th chebfun2 objcect that corresponds to
% a slice of F.
% length(F.slices{i}.rows) = the number of entries in each vector that
% corresponds to a row of the i-th slice of F.
% length(F.slices{i}.cols) = the number of entries in each vector that
% corresponds to a column of the i-th slice of F.
% length(F.pivotValues) = number of pivot values in the slice decomposition
% which is the same as F.rank.
% size(F.pivotLocations,1)*size(F.pivotLocations,2) = size of the matrix
% that contains the pivot locations of the slice decomposition of F.
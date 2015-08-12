function [f,h] = extractPole(f)
% EXTRACTPOLE  Removes from F the term accounting for the non-zero poles
%
% [G,H] = EXTRACTPOLE(F) if F is non-zero at the poles, then this function
% returns SPHEREFUNs G and H, such that H is rank 1, G = F - H, and G is
% zero at the poles. If F is zero at the poles then G is an empty
% SPHEREFUN.
%

if ~isempty(f) && f.nonZeroPoles
    h = f;
    % The col and row for the pole is always assumed to be at the first
    % index
    h.cols = f.cols(:,1);
    h.rows = f.rows(:,1);
    h.pivotValues = f.pivotValues(1);
    h.idxPlus = 1;
    h.idxMinus = [];
    h.pivotIndices = f.pivotIndices(1,:);
    h.pivotLocations = f.pivotLocations(1,:);
    
    % Now remove these rows and columns from f
    
    if size(f.cols,2) > 1
        f.cols = f.cols(:,2:end);
        f.rows = f.rows(:,2:end);
    else     % Bug in chebfun requires checking for these cases.
        f.cols = [];
        f.rows = [];
    end
    f.pivotValues = f.pivotValues(2:end);
    f.idxPlus = f.idxPlus(2:end)-1;
    f.idxMinus = f.idxMinus - 1;
    f.pivotIndices = f.pivotIndices(2:end,:);
    f.pivotLocations = f.pivotLocations(2:end,:);
    f.nonZeroPoles = 0;
else
    h = spherefun([]);
end

end

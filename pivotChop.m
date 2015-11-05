function f = pivotChop(f)
piv = chebtech2();
tempPivots = f.pivotValues;
if numel(tempPivots) < 17 % Increase the number of pivots to 17 if less.
    tempPivots(end+1:17) = tempPivots(end);
end
piv.coeffs = tempPivots';
piv = simplify(piv);
Rank = min(numel(piv.coeffs),numel(f.pivotValues)); % If the rank was already good enough make no changes.
f = chebfun2(f.cols(:,1:Rank)*diag(1./f.pivotValues(1:Rank))*f.rows(:,1:Rank)');
end

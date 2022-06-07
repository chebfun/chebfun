function [relTol, absTol] = getTol(M, pseudoLevel, tolOld,domDiff)
% getTol estimates the condition number of f based on the values stored in
% M, it provides updated tolerances for the constructor
% see https://github.com/chebfun/chebfun/issues/1491
relTol = 2*size(M,1)^(4/5) * pseudoLevel;
vscale = max(abs(M(:)));
cheb = @(i,n) -cos((i-1).*pi/(n-1));
points = 1:size(M,1);
points = cheb(points, size(M,1));
gradNorms = zeros([1,size(M,1)]);
for i = 1:size(M,2)
    gradNorms(i) = max(abs(diff(M(:,i)) ./ diff(points)'));
end
gradNorms = max(gradNorms);
absTol = max(max(domDiff.*gradNorms), vscale) * relTol;
absTol = max([absTol, tolOld, pseudoLevel]);
end
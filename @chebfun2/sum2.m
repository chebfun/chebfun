function out = sum2(F)

out = sum(F.cols)*diag(1./F.pivotValues)*sum(F.rows).';

end
function z = fevalm(F, x, y)

zCol = feval(F.cols, y(:));
zRow = feval(F.rows, x(:));

z = zCol*diag(1./F.pivotValues)*zRow.';

end
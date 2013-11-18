function z = feval(F, x, y)

sx = size(x);
sy = size(y);

if ( min(sx) > 1 && all(sx == sy) )
    if ( rank(x) == 1 && rank(y) == 1 )
        x = x(1,:);
        y = y(:,1);
    end
end

zCol = feval(F.cols, y(:));
zRow = feval(F.rows, x(:));

z = zCol*diag(1./F.pivotValues)*zRow.';

end
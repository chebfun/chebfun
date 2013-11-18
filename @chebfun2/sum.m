function f = sum(F, dim)

if ( nargin == 1 )
    dim = 1;
end

if ( dim == 1 )
    f = F.rows * (sum(F.cols)*diag(1./F.pivotValues)).';
    f = chebfun({bndfun(f)}).';
elseif ( dim == 2 )
    f = F.cols * (diag(1./F.pivotValues)*sum(F.rows).');
    f = chebfun({bndfun(f)});
else 
    error('chebfun3 doesn''t exist yet; fool.');
end


end
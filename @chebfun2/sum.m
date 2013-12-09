function f = sum( F, dim )

if ( nargin == 1 )
    dim = 1;
end

if ( dim == 1 )
    f = F.rows * (sum(F.cols) * diag(1./F.pivotValues)).';
elseif ( dim == 2 )
    f = F.cols * (diag(1./F.pivotValues) * sum(F.rows).');
else 
    error('CHEBFUN2:SUM:unknown', ...
          'Undefined function ''sum'' for that dimension');
end
end
function F = diff(F, k, dim)

if ( nargin < 2 )
    k = 1;
end
if ( nargin < 3 )
    dim = 1;
end

if ( dim == 1 )
    F.cols = diff(F.cols, k);
elseif ( dim == 2 )
    F.rows = diff(F.rows, k);
else 
    error('chebfun3 doesn''t exist yet; fool.');
end


end
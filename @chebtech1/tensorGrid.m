function out = tensorGrid(N, dom)

dim = numel(N);
if ( dim == 1 )
    m = N(1);
    out = chebpts(m, dom(1:2), 1);
elseif ( dim == 2 )
    m = N(1);
    n = N(2);
    out = cell(dim, 1);
    xx = chebpts(m, dom(1:2), 1);
    yy = chebpts(n, dom(3:4), 1);
    [xx, yy] = ndgrid(xx, yy);
    out{1} = xx;
    out{2} = yy;
else
    m = N(1);
    n = N(2);
    p = N(3);
    out = cell(dim, 1);
    xx = chebpts(m, dom(1:2), 1);
    yy = chebpts(n, dom(3:4), 1);
    zz = chebpts(p, dom(5:6), 1);
    [xx, yy, zz] = ndgrid(xx, yy, zz);
    out{1} = xx;
    out{2} = yy;
    out{3} = zz;
end

end
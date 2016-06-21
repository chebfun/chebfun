function out = tensorGrid(N, dom)

dim = numel(N);
if ( dim == 1 )
    m = N(1);
    out = trigpts(m, dom(1:2));
elseif ( dim == 2 )
    m = N(1);
    n = N(2);
    out = cell(dim, 1);
    xx = trigpts(m, dom(1:2));
    yy = trigpts(n, dom(3:4));
    [xx, yy] = ndgrid(xx, yy);
    out{1} = xx;
    out{2} = yy;
else
    m = N(1);
    n = N(2);
    p = N(3);
    out = cell(dim, 1);
    xx = trigpts(m, dom(1:2));
    yy = trigpts(n, dom(3:4));
    zz = trigpts(p, dom(5:6));
    [xx, yy, zz] = ndgrid(xx, yy, zz);
    out{1} = xx;
    out{2} = yy;
    out{3} = zz;
end

end
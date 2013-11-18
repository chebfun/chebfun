function varargout = contour(F)

r = length(F);
[m, n] = length(F);

x = chebpts(round(4*pi*m), F.domain([1 2]));
y = chebpts(round(4*pi*n), F.domain([3 4]));
[xx, yy] = meshgrid(x, y);
zz = fevalm(F, x, y);

h = contour(xx, yy, zz);

if ( nargout > 0 )
    varargout{1} = h;
end

end
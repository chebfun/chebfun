function varargout = plot(F)

r = length(F);
[m, n] = length(F);

m = 50;
n = 50;

x = chebpts(round(4*pi*m), F.domain([1 2]));
y = chebpts(round(4*pi*n), F.domain([3 4]));
[xx, yy] = meshgrid(x, y);
zz = fevalm(F, x, y);

h = surf(xx, yy, zz, 'facecolor','interp', 'edgealpha', 0);

if ( nargout > 0 )
    varargout{1} = h;
end

end
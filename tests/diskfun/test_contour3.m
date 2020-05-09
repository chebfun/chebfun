function pass = test_contour3( pref )
% Test contour3

if ( nargin == 0 )
    pref = chebfunpref;
end

f = diskfun(@(x,y) cos(cos(4*x).^2 + sin(5*y).^2));
x = -pi:.1:pi;
y = -1:.1:1;
[xx, yy] = meshgrid(x,y);

pass = 1;
try
   contour3(f)
   contour3(f, 5)
   contour3(f, [0.4 0.4])
   contour3(f, 'numpts', 100)
   contour3(f, 'pivots', 'r.-')
   contour3(xx, yy, f)
catch ME
    pass = 0;
end

close all

end

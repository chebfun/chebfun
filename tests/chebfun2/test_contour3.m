function pass = test_contour3( pref )
% Test contour3

if ( nargin == 0 )
    pref = chebfunpref;
end

f = chebfun2(@(x,y) cos(x.*y));
x = -1:.1:1;
[xx, yy] = meshgrid(x);

pass = 1;
try
   contour3(f)
   contour3(f, 5)
   contour3(f, [0.8 0.8])
   contour3(f, 'numpts', 100)
   contour3(f, 'pivots', 'r.-')
   contour3(xx, yy, f)
catch ME
    pass = 0;
end

close all

end

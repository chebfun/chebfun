function pass = test_contour3( pref )
% Test contour3

if ( nargin == 0 )
    pref = chebfunpref;
end

f = spherefun.sphharm(4, 3);

pass = 1;
try
   contour3(f)
   contour3(f, 5)
   contour3(f, [0.3 0.3])
   contour3(f, 'numpts', 100)
catch ME
    pass = 0;
end

close all

end

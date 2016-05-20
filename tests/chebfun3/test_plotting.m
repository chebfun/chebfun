function pass = test_plotting()
% Check that the very basic plotting commands do not crash.

f = chebfun3(@(x,y,z) exp(cos(10*x.*y.*z)));
hold off

% Non-GUI plots first
plot(f),                     j = ishold;
slice(f, 'noslider'),        j = j + ishold;
slice(f, 0.5, -0.3, 0.9),    j = j + ishold;
isosurface(f, 'noslider'),   j = j + ishold;
isosurface(f, [0.5, -0.6]),  j = j + ishold;
scan(f),                     j = j + ishold;
scan(f, 'hold'),             j = j + ishold;
scan(f, 1),                  j = j + ishold;
scan(f, 1, 'hold'),          j = j + ishold;
scan(f, 2),                  j = j + ishold;
scan(f, 2, 'hold'),          j = j + ishold;
scan(f, 3),                  j = j + ishold;
scan(f, 3, 'hold'),          j = j + ishold;
close all
if ( j == 0 )
    pass(1) = 1;
else
    pass(1) = 0;
end

% GUI plots now:
isosurface(f),               j = j + ishold;
slice(f),                    j = j + ishold;
surf(f),                     j = j + ishold;
close all
if ( j == 0 )
    pass(2) = 1;
else
    pass(2) = 0;
end

% rank-1 and off the default domain
f = chebfun3(@(x,y,z) x.*y.*z, [-1 2 -1 2 -1 2]);
plot(f)
surf(f)
slice(f)
isosurface(f)
close all
pass(3) = 1;

hh = @(x,y,z) .75*exp(-((9*x-2).^2 + (9*y-2).^2 + (9*z-2).^2)/9) + ...
    .75*exp(-((9*x+1).^2)/49 - (9*y+1) - (9*z+1)/10) + ...
    .5*exp(-((9*x-7).^2 + (9*y-3).^2)/4 + (9*z-7).^2) - ...
    .2*exp(-(9*x-4).^2 - (9*y-7).^2 - (9*z-7).^2);
h = chebfun3(hh);
plot(h)
close all
pass(4) = 1;

end
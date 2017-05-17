function pass = test_plotting( pref )
% Check that the very basic plotting commands do not crash in diskfunv

F = grad(diskfun(@(x,y) x.*cos(y)));
try
    hold off
    quiver(F),         j = ishold;
    close all
    if ( j == 0 )
        pass(1) = 1;
    else
        pass(1) = 0;
    end
catch
    pass(1) = 0;
end


end
function pass = test_plotting( pref )
% Check that the very basic plotting commands do not crash in spherefunv

F = grad(spherefun(@(x,y,z) x.*cos(y.*z)));
try
    hold off
    quiver(F),         j = ishold;
    quiver3(F),        j = ishold;
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
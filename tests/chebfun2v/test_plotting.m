function pass = test_plotting( pref )
% Check that the very basic plotting commands do not crash in Chebfun2v

F = chebfun2v(@(x,y) x,@(x,y) y);
G = chebfun2v(@(x,y) x,@(x,y) y,@(x,y)y);
try
    hold off
    quiver(F),         j = ishold;
    quiver3(F),        j = ishold;
    surf(G),           j = ishold;
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
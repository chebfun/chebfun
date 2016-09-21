function pass = test_quiver3()
% Check that the quiver3 plotting command does not crash in Chebfun3v

F = chebfun3v(@(x,y,z) x, @(x,y,z) y, @(x,y,z) z);
try
    hold off
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
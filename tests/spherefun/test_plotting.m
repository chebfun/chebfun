function pass = test_plotting( pref )
% Check that the very basic plotting commands do not crash. 

f = spherefun(@(lam,th) 1-exp(sin(lam).*cos(lam).*sin(th).^2));
hold off
plot(f),        j = ishold;
surf(f),        j = j + ishold;
contour(f),     j = j + ishold;
close all
if ( j == 0 )
    pass(1) = 1;
else
    pass(1) = 0;
end

% Check that axis limits are set properly.
plot(f)
ax = axis;
pass(2) = all( (ax-[-1 1 -1 1 -1 1]) == 0 );
close all

contour(f)
ax = axis;
pass(3) = all( (ax-[-1 1 -1 1 -1 1]) == 0 );
close all

% Simple testing that options don't crash
plot(f,'.-')
pass(4) = 1;
close all

contour(f,[0 0],'k-')
pass(5) = 1;
close all

end
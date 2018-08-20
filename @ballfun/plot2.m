function plot2(f)
% PLOT2 Plot a BALLFUN function on the ballfun and its slices

clf

% Plot f on the plane X-Y
subplot(2,2,2);
plot(extract_diskfun(f,'x','y'))
colorbar
xlabel('X')
ylabel('Y')

% Plot f on the plane X-Z
subplot(2,2,3);
plot(extract_diskfun(f,'x','z'));
colorbar
xlabel('X')
ylabel('Z')

% Plot f on the plane Y-Y
subplot(2,2,4);
plot(extract_diskfun(f,'y','z'));
colorbar
xlabel('Y')
ylabel('Z')

% Plot f
subplot(2,2,1);
plot(f, 'Grid')
title('f')

set(gcf,'PaperPositionMode','auto','PaperPosition',[0 0 15 10])
end

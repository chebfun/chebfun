function varargout = cheb2logo
%CHEB2LOGO   The unofficial Chebfun2 logo.
%   CHEB2LOGO plots the Chebfun2 logo.
%
%   H = CHEB2LOGO returns a figure handle to the logo.
%
% See also CHEBLOGO.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% TODO:  Comment this file.

LW = 'LineWidth';
FS = 'FontSize';

dom = [-1 0.957];

figure();
f = restrict(chebpoly(10), dom);
plot(f, LW, 5);
hold on

t = -cos(pi*(2:8)'/10)*0.99;  % Chebyshev extrema (tweaked).
y = 0*t;
h = text(t, y, num2cell(transpose('chebfun')), FS, 28, ...
    'hor', 'cen', 'vert', 'mid');

axis([-1.02 1.2 -1.5 3])
x = chebfun('x', dom);
plot(x + .2, f + 2, LW, 5)
hold on

f2 = f + 2;
x2 = x + .2;

xx = linspace(dom(1), dom(2), 1000).';
ff = [f(xx) f2(xx)];

close all
surf(repmat([0 1], length(xx), 1), -[xx xx], ff)
shading interp
colormap winter

fs = 18;

h1 = text(y, -t, num2cell(transpose('chebfun')), FS, fs, ...
    'hor', 'cen', 'vert', 'mid') ;
set(gca, 'view', [ -74.500000000000000 -32.000000000000000]);

t = linspace(t(1), t(end), length(t));
h2 = text(.5*(t + 1) - .05, -y - .98, (t + 1) - .23, ...
    num2cell(transpose('')), FS, fs, 'hor', 'cen', 'vert', 'mid', ...
    'Rotation', 10);

t = t(2:4);
h3 = text(.5*(t + 1) - .02, -y(1:3) - 1, (t + 1) + .028, ...
    num2cell(transpose('two')), FS, fs, 'hor', 'cen', 'vert', 'mid', ...
    'Rotation',10);

set(gca, 'view', [-72.5000 -20.0000]);
axis off
set(gcf, 'color', 'w');
set(gca, 'pos', [0 0 1 1]);

hold on
C = [0.2 0.2 1];
lw = 3;

plot3(0*xx, -xx, f(xx), LW, lw, 'Color', C)
plot3(0*xx + 1, -xx, f2(xx), 'b', LW, lw, 'Color', C)
plot3([0 1], -dom(2)*[1 1], f(xx(end)) + [0 2], 'b', LW, lw, 'Color', C)
plot3([0 1], [1 1], [1 3], 'b', LW, lw, 'Color', C)
colormap([1 1 1 ; .5 .5 .5]);

flist = listfonts;
k = strmatch('Rockwell', flist);  % 1st choice
if ( ~isempty(k) )
    set(h1, 'fontname', flist{k(1)})
    set(h2, 'fontname', flist{k(1)})
    set(h3, 'fontname', flist{k(1)})
end

set(gca, 'xlim', [-.1 1])
set(gcf, 'position', [440 480 560 240]);

oldscreenunits = get(gcf, 'Units');
oldpaperunits = get(gcf, 'PaperUnits');
oldpaperpos = get(gcf, 'PaperPosition');

set(gcf, 'Units', 'pixels');
scrpos = get(gcf, 'Position');

newpos = scrpos/100;
set(gcf, 'PaperUnits', 'inches', 'PaperPosition', newpos);
% print -dpng andbeyond_white
drawnow

set(gcf, 'Units', oldscreenunits, 'PaperUnits', oldpaperunits, ...
    'PaperPosition', oldpaperpos)

set(gcf, 'unit', 'cent', 'menuBar', 'none', ...
    'name', 'Chebfun2 logo', 'numbertitle', 'off', 'color', 'w')

if ( nargout == 1 )
    varargout = {h};
end

end

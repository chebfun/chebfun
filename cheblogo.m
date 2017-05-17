function varargout = cheblogo()
%CHEBLOGO   Plot the Chebfun logo.
%   CHEBLOGO plots the Chebfun logo.
%
%   F = CHEBLOGO returns a CHEBFUN of the Chebfun logo.

% Copyright 2017 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% Make a CHEBFUN of the logo:
f = chebpoly(10);
dom = [-1, .957];
f = restrict(f, dom);
x = chebfun('x', dom);

if ( nargout > 0 )
    % Export the logo:
    varargout{1} = f;
    return
end

figure
% Plot the shadow:
plot(x+.015, f-.075, 'color', .7*[1 1 1], 'LineWidth', 5);
hold on
% Plot the curve:
plot(f, 'b', 'LineWidth', 5)

% Plot the text:
t = - cos(pi*(2:8)'/10) *0.99;            % cheb extrema (tweaked)
y = 0*t; 
h = text( t, y, num2cell(transpose('chebfun')), ...
  'FontSize', 28, 'hor', 'cen', 'vert', 'mid') ;

% Choose a nice font:
flist = listfonts;
k = strmatch('Rockwell', flist);          % 1st choice
k = [k ; strmatch('Luxi Serif', flist)];  % 2nd choice
k = [k ; strmatch('Times', flist)];       % 3rd choice
if ( ~isempty(k) ) 
    set(h, 'FontName', flist{k(1)});
end

% Adjut the window size, etc.:
axis([-1.05 1 -1.8 1.8]), axis off
set(gca, 'pos', [0 0 1 1])
un = get(0, 'unit'); 
set(0, 'unit', 'cent')
ssize = get(0, 'screensize');  
set(0, 'unit', un)
set(gcf, 'papertype', 'A4', 'paperunit', 'cent', 'paperpos', [4.49 12.83 12 4])
pos = [ (ssize(3)-12)/2 (ssize(4)-4)/2 12 4];
set(gcf, 'unit', 'cent', 'pos', pos, 'menuBar', 'none', ...
    'name', 'Chebfun logo', 'numbertitle', 'off', 'color', 'w')

end

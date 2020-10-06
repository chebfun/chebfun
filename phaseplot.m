function phaseplot(f,ax)  % Plot phase portrait.
%PHASEPLOT   Phase (= argument) plot of a complex function
%   PHASEPLOT(F), where F is a function handle, draws a phase plot of F(Z) in
%   the complex plane.  As arg(f(z)) ranges over [0,2pi] the colors run
%   cyan -> blue -> red -> yellow -> green -> cyan.
%   
%   The existing axes are used, or [-1 1 -1 1] if there are no existing axes.
%   PHASEPLOT(F,[A B C D]) uses the axes [A B C D].
%
% Examples:
%   
%   phaseplot(@(z) z.^2), axis square
%   phaseplot(@(z) 1./z), axis square
%   phaseplot(@(z) besselj(6,z),[-12 12 -5 5]), axis equal
%
%   plot(chebfun('exp(1i*s)',[-pi,pi]),'k')
%   axis([-2 2 -2 2]), axis square
%   r = padeapprox(@(z) sqrt(1-z.^3),24,24);
%   phaseplot(r)
%
% Reference:
%
%   E. Wegert, Visual Complex Functions: An Introduction with Phase Portraits,
%   Springer Basel, 2012,
%
% See also CHEBFUN2/PLOT.

% Copyright 2020 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if nargin == 2
    figure
    axis(ax)
elseif isempty(get(gcf, 'CurrentAxes'))
    ax = [-1 1 -1 1];
    axis(ax)
else
    ax = axis;
end
holdstate = ishold;
hold on
x = linspace(ax(1), ax(2), 300);
y = linspace(ax(3), ax(4), 300);
[xx,yy] = meshgrid(x,y);
zz = xx+1i*yy;
h = surf(real(zz), imag(zz), -ones(size(zz)), angle(f(zz)));
set(h, 'EdgeColor','none');
caxis([-pi pi])
colormap hsv(600)
view(0,90)
if holdstate == 0
    hold off
end

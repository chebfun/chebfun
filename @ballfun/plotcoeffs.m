function varargout = plotcoeffs(f, varargin)
%PLOTCOEFFS   Display the coefficients of the columns, rows and tubes.
%   PLOTCOEFFS(F) plots the maximum coefficients in the column, row and tube
%   directions on a semilogy scale of the underlying basis. It returns three
%   figures one for the columns, one for the rows and one for the tubes.
%
% See also CHEBFUN/PLOTCOEFFS and CHEBFUN2/PLOTCOEFFS.

% Copyright 2019 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Empty check.
if ( isempty(f) )
    varargout = { [] }; 
    return
end

% Store the hold state of the current axis:
holdState = ishold;

% Get the tensor of coefficients of f
cfs = f.coeffs;

% Maximum Chebyshev coefficients along the radial direction
r_cfs = max(max( abs(cfs), [], 2), [], 3);
rr = chebfun(r_cfs, 'coeffs');
% Maximum Fourier coefficients along the azimuthal direction
l_cfs = max(max( abs(cfs), [], 1), [], 3);
l_cfs = l_cfs(:);
ll = chebfun(l_cfs, [-pi,pi], 'coeffs', 'trig');
% Maximum Fourier coefficients along the polar direction
t_cfs = max(max( abs(cfs), [], 1), [], 2);
t_cfs = t_cfs(:);
tt = chebfun(t_cfs, [-pi,pi], 'coeffs', 'trig');

% PLOTCOEFFS of cols:
ax1 = subplot(1, 3, 1);
plotcoeffs(rr, varargin{:}); 
ylim1 = ylim(gca);
% Remove labels from 1D plotcoeff: 
ylabel(gca, 'Magnitude of coefficient') 
title('Cols (r)')  

% PLOTCOEFFS of rows:
ax2 = subplot(1, 3, 2);
plotcoeffs(ll, varargin{:}); 
ylim2 = ylim(gca);
ylabel(' ')
title('Rows (lambda)')

% PLOTCOEFFS of tubes:
ax3 = subplot(1, 3, 3);
plotcoeffs(tt, varargin{:}); 
ylim3 = ylim(gca);
ylabel(' ')
title('Tubes (theta)')

% Find a proper ylim for all the three subplots:
yLims = [ylim1; ylim2; ylim3];
ylimNew = [min(yLims(:, 1)), max(yLims(:, 2))];

% Set the ylim of the plots again to be similar.
ylim(ax1, ylimNew);
ylim(ax2, ylimNew);
ylim(ax3, ylimNew);

% Return hold state to what it was before:
if ( ~holdState )
    hold off
end
end



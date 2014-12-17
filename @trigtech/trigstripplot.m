function varargout = trigstripplot(u, varargin)
%TRIGSTRIPPLOT   Plot the strip of analyticity.
%   TRIGSTRIPPLOT(U) plots estimated strip of analyticity in the complex plane
%   for U. The width of the strip is 2*a(k)=1/N(k)/pi*log(4/EPS+1), where EPS is
%   the EPSLEVEL of U and N(k) is the number of Fourier modes.
%
%   TRIGSTRIPPLOT(U, EPS) allows a user-specified EPS.
%
%   TRIGSTRIPPLOT(U, ..., S) allows plotting options to be passed. For
%   example, for black lines one may write TRIGSTRIPPLOT(U, 'k-').
%
%   H = TRIGSTRIPPLOT(U) returns a handle H to the figure.
%
%   Example:
%       u = trigtech(@(x) 1./(1 + sin(pi*x).^2));
%       trigstripplot(u, 'r--');

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( size(u, 2) > 1 )
    error('CHEBFUN:TRIGTECH:trigstripplot:quasi', ...
        ['TRIGSTRIPPLOT does not support array-valued TRIGTECH objects.']);
end

if ( isempty(u) )
    h = plot([]);
    % Output the axis handle.
    if ( nargout ~= 0 )
        varargout = {h};
    end
    return
end

% Parse the inputs.
[ee, args] = parseInputs(varargin{:});
if ( isnan(ee) )
    ee = u.epslevel;
end

% Compute the intercept of the strip on the imaginary axis
M = length(u);
if ( mod(M,2) == 1 )
    N = (M-1)/2;
else
    N = M/2-1;
end
a = 1i/N/pi*log(4/ee+1);

UK = [-1 + a ; 1 + a ; nan ; -1 - a ; 1 - a];

holdState = ishold();

% Plot the ellipses.
h = plot(UK, args{:}); hold on
clr = get(h, 'Color');

% Plot the interval (with ticks).
dom = [-1 1];
h2 = plot(dom, 0*dom, args{:});
set(h2, 'color', [0 0 0], 'marker', '+');

UK  = [-1+a; 1+a;1-a; -1-a;-1+a];
h3 = fill(real(UK), imag(UK), clr, 'EdgeColor', 'none', 'FaceAlpha', 0.2);

h = [h; h2; h3];

if ( ~holdState )
    hold off
end

% Output the axis handle.
if ( nargout ~= 0 )
    varargout = {h};
end

end

function [ee, args] = parseInputs(varargin)

% Default options
ee = NaN;               % Default EPS
args = {};              % Additional plotting args.

% Sort out the inputs.
if ( nargin >= 1 )
    % Check arguments for EPS and K, if they exist.
    if ( isnumeric(varargin{1}) )
        ee = varargin{1};
        varargin(1) = [];
    end
    if ( (numel(varargin) >= 1) )
        args = varargin;
    end
end

end
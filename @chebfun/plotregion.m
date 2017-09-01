function varargout = plotregion(u, varargin)
%PLOTREGION   Plot estimated regions of analyticity.
%   PLOTREGION(U) plots the estimated region of analyticity in the complex
%   plane for each piecewise part of U. If U is a Chebyshev-based chebfun, this 
%   is the ellipse with foci at the endpoints of U.domain and semi-minor
%   and semi-major axes summing to rho(k) = C*exp(abs(log(EPS))/N(k)),
%   where C is the appropriate scaling for the interval
%   [U.ends(k) U.ends(k+1)] and EPS is the machine precision. If U
%   is a trigonometric-based periodic chebfun, it is the strip symmetric
%   about the real axis with half-width log(1/EPS)/(pi*N).
%
%   PLOTREGION(U, EPS) allows a user-specified EPS.
%
%   PLOTREGION(U, K) and PLOTREGION(U, EPS, K) plot the estimated regions of 
%   analyticity for the funs of U indexed by the vector K.
%
%   PLOTREGION(U, ..., S) allows plotting options to be passed. For example, for
%   black lines one may write PLOTREGION(U, 'k-').
%
%   PLOTREGION(U, ..., 'legends', 0) prevents the legends being displayed on the
%   plot.
%
%   PLOTREGION(U, ..., 'numpts', N) plots each ellipse using N points without 
%   affecting the strips if any.
%
%   H = PLOTREGION(U) returns a handle H to the figure.
%
%   Example:
%       u = chebfun({@sin, @cos, @tan, @cot}, [-2, -1, 0, 1, 2]);
%       plotregion(u)

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( numColumns(u) > 1 )
    error('CHEBFUN:CHEBFUN:plotregion:quasi', ...
        ['CHEBELLPISEPLOT does not support array-valued CHEBFUN objects ' ...
         'or quasimatries.']);
end

if ( isempty(u) )
    h = plot([]);
    % Output the axis handle.
    if ( nargout ~= 0 )
        varargout = {h};
    end
    return
end

% Simplify first:
u = simplify(u);

% Parse the inputs.
[k, ee, numpts, legends, args] = parseInputs(varargin{:});
if ( isnan(ee) )
    ee = eps;
end
[lineStyle, pointStyle, jumpStyle, deltaStyle, args] = ...
    chebfun.parsePlotStyle(args{:});

% If k==0 (as by default), we will plot all funs.
if ( k == 0 )
    k = 1:length(u.funs);
end

% Transpose U properly.
if ( u.isTransposed )
    u = u';
end

% Error if index exceeds dimensions.
if ( any(k > length(u.funs)) )
    error( 'CHEBFUN:CHEBFUN:plotregion:outOfBounds', ...
        'Input chebfun has only %d pieces', length(u.funs) );
end

% A cell array of ellipse coordinates to be plotted.
UK = {};
AUX = {};
xl = [];
yl = [];
for j = k
    uk = u.funs{j};
    endsk = uk.domain;
    if ( any(isinf(endsk)) )
        error( 'CHEBFUN:CHEBFUN:plotregion:unboundedDomain', ...
            ['Plot of analyticity region is not supported for function on' ...
            'unbounded domain.']);
    end
    data = plotregionData(uk, ee, numpts);
    ek = .5*sum(endsk) + .5*diff(endsk)*data.boundary;
    UK = [UK, {real(ek), imag(ek)}, args{:}]; % Add the variable args.
    au = .5*sum(endsk) + .5*diff(endsk)*data.auxiliary;
    AUX = [AUX, {real(au), imag(au)}];
    xl = [xl .5*sum(endsk) + .5*diff(endsk)*data.xlim];
    yl = [yl .5*diff(endsk)*data.ylim];
end

xl = [min(xl) max(xl)];
yl = [min(yl) max(yl)];

holdState = ishold();

% Plot the ellipses.
h = plot(UK{:}, lineStyle{:});
cl = get(h, {'Color'});
hold on

% Plot the auxiliary lines:
for j = 1:numel(cl)
    hh = plot(AUX{2*(j-1)+1}, AUX{2*(j-1)+2}, ':');
    set(hh, 'Color', cl{j})
end

% Add the legend.
if ( legends ) && ( j > 1 )
    legend(int2str(k.'))
end

% Plot the interval (with ticks).
dom = u.domain;
h2 = plot(dom, 0*dom, args{:}, lineStyle{:}, pointStyle{:});
set(h2, 'color', [0 0 0], 'marker', '+', 'LineStyle', '-');
h = [h ; h2];

xlim(xl)
ylim(yl)

if ( ~holdState )
    hold off
end

% Output the axis handle.
if ( nargout ~= 0 )
    varargout = {h};
end

end

function [k, ee, numpts, legends, args] = parseInputs(varargin)

% Default options
k = 0;                  % plot all funs by default
ee = NaN;               % Default EPS
numpts = 101;           % Number of points in plots
legends = 1;            % Display legends?             

isPosInt = @(v) abs(round(v)) == v;

% Sort out the inputs.
if ( nargin >= 1 )
    % Check arguments for EPS and K, if they exist.
    v1 = varargin{1};
    if ( isnumeric(v1) )
        if ( v1 == 1 )
            k = v1;
        else
            if ( isPosInt(v1) && ~isnumeric(varargin{2}) )
                k = v1;
            elseif ( isnumeric(varargin{2}) )
                ee = v1;
                k = varargin{2};
                varargin(2) = [];
            end
        end
        varargin(1) = [];
    end
    
    % Check named inputs, if they exist.
    j = 1;
    while ( j < length(varargin) )
        if ( strcmpi(varargin{j}, 'numpts') )
            numpts = varargin{j+1}; 
            varargin(j:j+1) = [];
        elseif ( strcmpi(varargin{j}, 'legends') )
            legends = varargin{j+1}; 
            varargin(j:j+1) = [];
        else
            j = j+1;
        end
    end
end

args = varargin; % Additional plotting args.

end

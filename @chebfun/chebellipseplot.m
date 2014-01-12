function varargout = chebellipseplot(u, varargin)
%CHEBELLIPSEPLOT   Plot the Bernstein (aka Chebyshev) ellipses.
%   CHEBELLIPSEPLOT(U) plots Bernstein ellipses in the complex plane for each
%   piecewise part of U, with foci at points in U.domain and semi-minor and
%   major axes summing to rho(k) = C*exp(abs(log(EPS))/N(k)), where C is the
%   appropriate scaling for the interval [U.ends(k) U.ends(k+1)] and EPS is the
%   EPSLEVEL of U.
%
%   CHEBELLIPSEPLOT(U, EPS) allows a user-specified EPS.
%
%   CHEBELLIPSEPLOT(U, K) and CHEBELLIPSEPLOT(U, EPS, K) plot ellipses for
%   the funs of U indexed by the vector K. If U is a quasimatrix, only
%   the first column/row is considered.
%
%   CHEBELLIPSEPLOT(U, ..., S) allows plotting options to be passed. For
%   example, for black lines one may write CHEBELLIPSEPLOT(U, 'k-').
%
%   CHEBELLIPSEPLOT(U, ..., 'legends', 0) prevents the legends being
%   displayed on the plot.
%
%   CHEBELLIPSEPLOT(U, ..., 'numpts', N) plots each ellipse using N points.
%
%   H = CHEBELLIPSEPLOT(U) returns a handle H to the figure.
%
%   Example:
%       u = chebfun({@sin, @cos, @tan, @cot}, [-2, -1, 0, 1, 2]);
%       chebellipseplot(u, sqrt(eps), '--');

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.

%%

if ( min(size(u)) > 1 )
    % TODO: Replace with "NUMCOLS(U) > 1" from downstream in feature-quasimatrix.
    error('CHEBFUN:chebellipseplot:quasi', ['CHEBELLPISEPLOT does not ', ... 
        'support array-valued CHEBFUN objects or  quasimatries.']);
end

if ( isempty(u) )
    h = plot([]);
    return
end

% Parse the inputs.
[k, ee, numpts, legends] = parseInputs(varargin{:});
if ( isnan(ee) )
    ee = epslevel(u);
end

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
    error( 'CHEBFUN:chebellipseplot:outOfBounds', ...
        'Input chebfun has only %d pieces', length(u.funs) );
end

% The unit circle.
c = exp(2*pi*1i*linspace(0, 1, numpts));

% A cell array of ellipse coordinates to be plotted.
UK = {};
for j = k
    uk = u.funs{j};
    endsk = uk.domain;
    rhok = exp(abs(log(ee)) / length(uk));
    ek = .5*sum(endsk) + .25*diff(endsk)*(rhok*c + 1./(rhok*c));
    UK = [UK, {real(ek), imag(ek)}, varargin{:}]; % Add the variable args.
end

holdState = ishold();
hold on

% Plot the ellipses.
h = plot(UK{:});

% Add the legend.
if ( legends ) && ( j > 1 )
    legend(int2str(k.'))
end

% Plot the interval (with ticks).
dom = u.domain;
h2 = plot(dom, 0*dom, varargin{:});
set(h2, 'color', [0 0 0], 'marker', '+');
h = [h ; h2];

if ( ~holdState )
    hold off
end

% Output the axis handle.
if ( nargout ~= 0 )
    varargout = {h};
end

end

function [k, ee, numpts, legends] = parseInputs(varargin)

% Default options
k = 0;                  % plot all funs by default
ee = NaN;               % Default EPS
numpts = 101;           % Number of points in plots
legends = 1;            % Display legends?

% Sort out the inputs.
if ( nargin > 1 )
    % Check arguments for EPS and K, if they exist.
    if ( isnumeric(varargin{1}) )
        if ( varargin{1} >= 1 )
            k = varargin{1};
        else
            ee = varargin{1};
            if ( (numel(varargin) > 1) && (isnumeric(varargin{2})) )
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

end

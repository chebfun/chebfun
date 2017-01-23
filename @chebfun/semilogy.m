function varargout = semilogy(varargin)
%SEMILOGY   Semi-log scale plot of a CHEBFUN.
%   SEMILOGY(...) is the same as PLOT(...), except a logarithmic (base 10) scale
%   is used for the Y-axis.
%
% See also PLOT, SEMILOGX, LOGLOG.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Standard CHEBFUN/PLOT():
[h{1:3}] = plot(varargin{:});

% Find the CHEBFUNS:
isCheb = cellfun(@(v) isa(v, 'chebfun'), varargin);
% Choose a TOL based on vscale:
tol = max(cellfun(@(f) max(eps*vscale(f)), varargin(isCheb)));
% Loop over the different plot components:
for j = 1:numel(h)
    % Get the y data:
    yData = get(h{j}, 'yData');
    % Ensure it's a cell:
    if ( ~iscell(yData) )
        yData = {yData};
    end
    % Loop over each cell:
    for k = 1:numel(yData)
        if ( min(yData{k}(:)) > -tol )
            % If it's only a little negative, then take the absolute value:
            yData{k} = abs(yData{k});
        end
    end
    for k = 1:numel(h{j})
        % Put it back in h
        set(h{j}(k), 'yData', yData{k});
    end
end

% Set the YScale to be logarithmic:
set(gca, 'YScale', 'log');     

% Output handle if requested:
if ( nargout > 0 )
    varargout = {h};
end

end

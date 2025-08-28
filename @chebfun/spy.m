function spy(A, varargin)
%SPY   Visualise a CHEBFUN.
%   SPY(F) plot a visualisation of the CHEBFUN F. Break points are shown as gaps
%   in the line plot.
%
%   SPY(F, 'LineSpec') uses the color and marker from the line specification
%   string 'LineSpec' (See PLOT for possibilities).

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

hold off

% Grab the domain:
[a, b] = domain(A);
if ( isinf(a) )
    a = -1e16;
end
if ( isinf(b) )
    b = 1e16;
end
scl = (b-a)/100;
ee = [-scl ; NaN ; scl];

isTransposed = A(1).isTransposed;
numCols = numColumns(A);

% Loop over the columns:
for j = 1:numCols
    Aj = extractColumns(A, j);
    domj = Aj.domain;
    ss = repmat(domj, 3, 1) + repmat(ee, 1, length(domj));
    ss = ss(:);
    ss([2:3 end-1:end]) = [];
    ss([1, end]) = [a, b];
    jj = repmat(j, length(ss), 1);
    if ( isTransposed )
        [jj, ss] = deal(ss, jj);
    end
    if ( nargin > 1 ) 
        linespec = varargin{:};
    else
        linespec = ''; 
    end
    marks = {'.-', 'o', 'x', '+', '*', 's', 'd', 'v', '^', '<', '>', 'p', 'h'}; 
    reg = regexp(linespec, marks );
    matlab_blue = [0 .45 .74];
    if ( ~isempty( cell2mat( reg ) ) )
        % This allows the zero structure of a quasimatrix to be plotted
        rts = roots( Aj );
        rj = repmat(j, length(rts), 1);
        plot(jj, ss, varargin{:}, 'markersize', eps, 'color', matlab_blue)
        hold on 
        plot(rj, rts, varargin{:}, 'linestyle', 'none')
    else
        if ( nargin > 1 )
            plot(jj, ss, varargin{:}); 
        else
            plot(jj, ss, 'color', matlab_blue); 
        end
        hold on
    end

end

% Tidy the axes, etc.:
if ( ~isTransposed )  
    % Switch the plot direction:
    set(gca, 'ytick', [a, b], 'ydir', 'reverse')
    % Add some ticks:
    if ( numCols < 10 )
        set(gca, 'xtick', 1:numCols)
    else
        set(gca, 'xtick', [1 numCols]),
    end
    % Adjust the plot ratio:
    axis([0 numCols+1 a b])
    ar = get(gca, 'PlotBoxAspectRatio');
    ar(1) = .5*ar(1);
    set(gca, 'PlotBoxAspectRatio', ar)
else
    % Switch the plot direction:
    set(gca, 'xtick', [a b], 'ydir', 'reverse')
    % Add some ticks:
    if ( numCols < 10 )
        set(gca, 'ytick', 1:numCols)
    else
        set(gca, 'ytick', [1 numCols])
    end
    % Adjust the plot ratio:
    axis([a b 0 numCols+1])
    ar = get(gca, 'PlotBoxAspectRatio');
    ar(2) = .4*ar(2);
    set(gca, 'PlotBoxAspectRatio', ar)
end

hold off

end

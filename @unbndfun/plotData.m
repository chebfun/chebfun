function data = plotData(f, g)
%PLOTDATA    Useful data values for plotting an UNBNDFUN object.
%   DATA = PLOTDATA(F) returns a cell array of data values that can be used for
%   plotting F. In particular, DATA is a 4x1 cell array of the form {xLine,
%   fLine, xPoints, fPoints}, where xLine-fLine are a data pair for plotting the
%   continuous function F and xPoints-fPoints are the data pair for plotting
%   values of F on the underlying Chebyshev grid.
%
% See also PLOT.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( nargin == 1 || isempty(g) )
    % Get the data from the ONEFUN:
    data = plotData(f.onefun);
    
    % Map the 'x' data using f.mapping.for:
    data.xLine = f.mapping.for(data.xLine);
    data.xPoints = f.mapping.for(data.xPoints);
    
    %% Figure out the xlim:
    data.xLim = get(f, 'domain');
    
    % Size of the window:
    window = 10;
    
    % [TODO]: We need to find a better way to determine the center and the
    % width of the window.
    
    % Center the window at the origin.
    center = 0;
    
    % If the left endpoint is -Inf:
    if ( isinf(data.xLim(1)) )
        data.xLim(1) = center - window;
    end
    
    % If the right endpoint is Inf:
    if ( isinf(data.xLim(2)) )
        data.xLim(2) = center + window;
    end
    
    % Sort out the jumps:
    data.xJumps = [f.domain(1) ; NaN ; f.domain(2)];
    data.yJumps = getJumps(f, data.yLine);
    
end

end

function jumps = getJumps(f, fLine)
    lvalF = get(f, 'lval');
    rvalF = get(f, 'rval');
    
    % Deal with functions which blow up:
    ind = isinf(lvalF);
    lvalF(ind) = fLine(2, ind);
    ind = isinf(rvalF);
    rvalF(ind) = fLine(end-1, ind);
    
    myNaN = nan(size(lvalF));
    jumps = [lvalF ; myNaN ; rvalF];
end
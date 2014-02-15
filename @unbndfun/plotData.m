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
    
    if ( isinf(f) )
        
        % If F is infinite, figure out where the mininum of ABS(F) takes place:
        [ignored, idx] = min(abs(get(f, 'values')));
        pts = get(f, 'points');
        minLoc = f.mapping.for(pts(idx));
        
        % If the left endpoint is -Inf:
        if ( isinf(data.xLim(1)) )
            if ( isfinite(minLoc) )
                data.xLim(1) = minLoc - window;
            else
                data.xLim(1) = data.xLim(2) - window;
            end
        end
        
        % If the right endpoint is Inf:
        if ( isinf(data.xLim(2)) )
            if ( isfinite(minLoc) )
                data.xLim(2) = minLoc + window;
            else
                data.xLim(2) = data.xLim(1) + window;
            end
        end
        
    else
        % If F is finite, figure out where the maxinum of ABS(F) takes place:
        [ignored, idx] = max(abs(get(f, 'values')));
        pts = get(f, 'points');
        maxLoc = f.mapping.for(pts(idx));
        
        % If the left endpoint is -Inf:
        if ( isinf(data.xLim(1)) )
            data.xLim(1) = maxLoc - window;
        end
        
        % If the right endpoint is Inf:
        if ( isinf(data.xLim(2)) )
            data.xLim(2) = maxLoc + window;
        end
        
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
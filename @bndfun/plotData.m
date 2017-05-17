function data = plotData(f, g, h)
%PLOTDATA   Useful data values for plotting a BNDFUN object.
%   DATA = PLOTDATA(F) returns a struct containing data that can be used for
%   plotting F. In particular, DATA.xLine and DATA.yLine are for plotting smooth
%   curves (usually passed to plot with '-') and DATA.xPoints and DATA.yPoints
%   contain the (x, F(x)) data used to represent F.
%
%   DATA = PLOTDATA(F, G) returns data for PLOT(F, G), i.e., (F(x), G(x)), and
%   DATA = PLOTDATA(F, G, H) returns data for plots of the form PLOT3(F, G, H).
%   In the latter case, DATA also contains the fields zLine and zPoints.
%
% See also PLOT.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Get the data from the ONEFUN:
if ( (nargin == 1) || isempty(g) )
    % PLOT(F):
    data = plotData(f.onefun);
    % Map the 'x' data using f.mapping.For:
    data.xLine = f.mapping.For(data.xLine);
    data.xPoints = f.mapping.For(data.xPoints);
    
    % Sort out the jumps:
    data.xJumps = [f.domain(1) ; NaN ; f.domain(2)];
    data.yJumps = getJumps(f, data.yLine);
    
    % Overwrite the current xLim if it is provided by plotData@singfun:
    data.xLim = f.domain;
    
elseif ( nargin == 2 )
    % PLOT(F, G):
    data = plotData(f.onefun, g.onefun);
    
    % Sort out the jumps:
    data.xJumps = getJumps(f, data.xLine);
    data.yJumps = getJumps(g, data.yLine);
    
else
    % PLOT(F, G, H):
    data = plotData(f.onefun, g.onefun, h.onefun);
    
    % Sort out the jumps:
    data.xJumps = getJumps(f, data.xLine);
    data.yJumps = getJumps(g, data.yLine);
    data.zJumps = getJumps(h, data.zLine);
    
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

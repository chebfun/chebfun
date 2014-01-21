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

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Get the data from the ONEFUN:
if ( nargin == 1 || isempty(g) )
    % PLOT(F):
    
    data = plotData(f.onefun);
    % Map the 'x' data using f.mapping.for:
    data.xLine = f.mapping.for(data.xLine);
    data.xPoints = f.mapping.for(data.xPoints);
   
    lval = get(f, 'lval');
    rval = get(f, 'rval');
    
    % Consider the jump values for an infinite FUN:
    data.xJumps = [f.domain(1) ; NaN ; f.domain(2)];
    
    ind = isinf(lval);
    if ( any( ind ) )
        lval(ind) = data.yLine(1, ind);
    end
    
    ind = isinf(rval);
    if ( ind )
        rval = data.yLine(end, ind);
    end
    
    myNaN = nan(size(lval));
    data.yJumps = [lval ; myNaN ; rval];
    
elseif ( nargin == 2 )
    % PLOT(F, G):
    
    data = plotData(f.onefun, g.onefun);
    
    lvalG = get(g, 'lval');
    rvalG = get(g, 'rval');
    
    % Consider the jump values for an infinite FUN:
    lvalF = get(f, 'lval');
    rvalF = get(f, 'rval');

    myNaN = nan(size(lvalF));
    data.xJumps = [lvalF ; myNaN ; rvalF];
    
    ind = isinf(lvalG);
    if ( ind )
        lvalG(ind) = data.yLine(1, ind);
    end
    
    ind = isinf(rvalG);
    if ( ind )
        rvalG(ind) = data.yLine(end, ind);
    end
    
    myNaN = nan(size(lvalG));
    data.yJumps = [lvalG ; myNaN ; rvalG];
    
else
    % PLOT(F, G, H):
    data = plotData(f.onefun, g.onefun, h.onefun);
    
end

end

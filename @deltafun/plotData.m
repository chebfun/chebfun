function data = plotData(f)
%PLOTDATA   Useful data values for plotting a SINGFUN object.
%   DATA = PLOTDATA(F) extracts PLOTDATA of the smooth part of F
%   and then scales it by the singular factors given in the EXPONENTS of F
%
% See also SMOOTHFUN/PLOTDATA.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

%%
% Get plot data from the smooth parth
data = plotData( f.smoothPart );
% Update extrapolated y-data 
x = data.xLine;
y = data.fLine;
y = y.*(x + 1).^f.exponents(1);
y = y.*(1 - x).^f.exponents(2);
data.fLine = y;

% Update sample point y-data
y = data.fPoints;
x = data.xPoints;
y = y.*(x + 1).^f.exponents(1);
y = y.*(1 - x).^f.exponents(2);
data.fPoints = y;

end



function data = plotData(f, g, h)
%PLOTDATA   Useful data values for plotting a CHEBFUN object.
%   OUT = PLOTDATA(F) returns a struct containing data that can be used for
%   plotting F. In particular, DATA.xLine and DATA.yLine are for plotting smooth
%   curves (usually passed to plot with '-'), DATA.xPoints and DATA.yPoints
%   contain the (x, F(x)) data used to represent F, and DATA.xJumps and
%   DATA.yJumps are the linear connections between discontinuous pieces.
%
%   OUT = PLOTDATA(F, G) is similar, but for plots of the form PLOT(F, G). (Note
%   that F and G are assumed to be real-valued CHEBFUN objects). Here OUT.xLine,
%   OUT.xPoints, and OUT.xJumps contain the data relating to F, and OUT.yLine,
%   OUT.yPoints, OUT.yJumps the data relating to G.
%
%   OUT(F, G, H) returns data for plots of the form PLOT3(F, G, H), where F, G,
%   and H are real-valued CHEBFUN objects. In this case, the OUT also contains
%   the fields zLine, zPoints, and zJumps, which contain the plotting data
%   relating to H
%
% See also PLOT, PLOT3.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Initialise the output structure:
data = struct('xLine', [], 'yLine', [], 'xPoints', [], 'yPoints', [], ...
    'xJumps', [], 'yJumps', []);

if ( nargin == 1 || isempty(g) )
    % PLOT(F)

    % Loop over each FUN for Line and Points data:
    for k = 1:numel(f.funs)
        % Get the data from the FUN:
        dataNew = plotData(f.funs{k});
        myNaN = NaN(1, size(dataNew.yLine, 2)); % Array of NaNs.
        % Insert a NaN (or array of NaNs) and append new data to array:
        data.xLine = [data.xLine ; NaN ; dataNew.xLine];
        data.yLine = [data.yLine ; myNaN ; dataNew.yLine];
        data.xPoints = [data.xPoints ; NaN ; dataNew.xPoints];
        data.yPoints = [data.yPoints ; myNaN ; dataNew.yPoints];
    end

    % Return NaNs if there are no jumps:
    data.xJumps = NaN;
    data.yJumps = myNaN;
    % Loop over each FUN for Jumps data:
    for k = 1:(numel(f.funs) - 1)
        data.xJumps = [data.xJumps ; NaN ; f.funs{k}.domain(2) ; ...
            f.funs{k+1}.domain(1)];
        data.yJumps = [data.yJumps ; myNaN ; get(f.funs{k}, 'rval') ; ...
            get(f.funs{k+1}, 'lval')];
    end
    
elseif ( nargin == 2 )
    % PLOT(F, G)

    [f, g] = overlap(f, g);
    
    % Loop over each FUN for Line and Points data:
    for k = 1:numel(f.funs)
        % Get the data from the FUN objects:
        dataNew = plotData(f.funs{k}, g.funs{k});
        myNaN = NaN(1, size(dataNew.yLine, 2)); % Array of NaNs.
        % Insert a NaN (or array of NaNs) and append new data to array:
        data.xLine = [data.xLine ; myNaN ; dataNew.xLine];
        data.yLine = [data.yLine ; myNaN ; dataNew.yLine];
        data.xPoints = [data.xPoints ; myNaN ; dataNew.xPoints];
        data.yPoints = [data.yPoints ; myNaN ; dataNew.yPoints];
    end

    % Loop over each FUN for Jumps data:
    for k = 1:(numel(f.funs) - 1)
        % Append [oldData, NaN, rval_k, lval_{k+1}]:
        data.xJumps = [data.xJumps ; myNaN ; get(f.funs{k}, 'rval') ; ...
            get(f.funs{k+1}, 'lval')];
        data.yJumps = [data.yJumps ; myNaN ; get(g.funs{k}, 'rval') ; ...
            get(g.funs{k+1}, 'lval')];        
    end

else
    % PLOT(F, G, H)
    
    % Initialise z storage:
    data.zLine = [];
    data.zPoints = [];
    data.zJumps = [];

    % [TODO]: Fix this once OVERLAP() is implemented.
    if ( any( f.domain ~= g.domain ) )
        [f, g] = overlap(f, g);
    end
    if ( any( g.domain ~= h.domain ) )
        [g, h] = overlap(g, h);
        [h, f] = overlap(h, f);
    end
    
    % Loop over each FUN for Line and Points data:
    for k = 1:numel(f.funs)
        % Get the data from the FUN objects:
        dataNew = plotData(f.funs{k}, g.funs{k}, h.funs{k});
        myNaN = NaN(1, size(dataNew.yLine, 2)); % Array of NaNs.
        % Insert a NaN (or array of NaNs) and append new data to array:
        data.xLine = [data.xLine ; myNaN ; dataNew.xLine];
        data.yLine = [data.yLine ; myNaN ; dataNew.yLine];
        data.zLine = [data.zLine ; myNaN ; dataNew.zLine];        
        data.xPoints = [data.xPoints ; myNaN ; dataNew.xPoints];
        data.yPoints = [data.yPoints ; myNaN ; dataNew.yPoints];
        data.zPoints = [data.zPoints ; myNaN ; dataNew.zPoints];        
    end

    
    % Loop over each FUN for Jumps data:
    for k = 1:(numel(f.funs) - 1)
        % Append [oldData, NaN, rval_k, lval_{k+1}]:
        data.xJumps = [data.xJumps ; myNaN ; get(f.funs{k}, 'rval') ; ...
            get(f.funs{k+1}, 'lval')];
        data.yJumps = [data.yJumps ; myNaN ; get(g.funs{k}, 'rval') ; ...
            get(g.funs{k+1}, 'lval')]; 
        data.zJumps = [data.zJumps ; myNaN ; get(h.funs{k}, 'rval') ; ...
            get(h.funs{k+1}, 'lval')];         
    end
    
end
    
    
end


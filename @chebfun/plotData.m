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

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Initialise the output structure:
data = struct('xLine', [], 'yLine', [], 'xPoints', [], 'yPoints', [], ...
    'xJumps', [], 'yJumps', [], 'xDeltas', [], 'yDeltas', [], ...
    'yDeltaBase', [], 'xLim', [Inf -Inf], 'yLim', [Inf -Inf], ...
    'defaultXLim', 1, 'defaultYLim', 1);

if ( nargin == 1 )
    % PLOT(F)

    % Loop over each FUN for Line and Points data:
    nFuns = numel(f.funs);
    
    for k = 1:nFuns
        % Get the data from the FUN:
        dataNew = plotData(f.funs{k});

        if ( ~any(ismember(fields(dataNew), 'xDeltas' )) )
            dataNew.xDeltas = NaN;
            dataNew.yDeltas = NaN;
            dataNew.yDeltaBase = NaN;
        end

        data.xLim = [min(dataNew.xLim(1), data.xLim(1)), ...
            max(dataNew.xLim(2), data.xLim(2))];
        data.yLim = [min(dataNew.yLim(1), data.yLim(1)), ...
            max(dataNew.yLim(2), data.yLim(2))];
        
        if ( k == 1 )
            dataNew.xJumps(1) = [];
            dataNew.yJumps(1,:) = [];
        end
        
        if ( k == nFuns )
            dataNew.xJumps(end) = [];
            dataNew.yJumps(end,:) = [];
        end
        
        if ( any(~isreal(dataNew.yLine)) )
            % Deal with complex-valued functions:
            
            % Assign x to be the real part, and y to be the imaginary part:
            dataNew.xLine = real(dataNew.yLine);
            dataNew.yLine = imag(dataNew.yLine);
            dataNew.xPoints = real(dataNew.yPoints);
            dataNew.yPoints = imag(dataNew.yPoints);
            dataNew.xJumps = real(dataNew.yJumps);
            dataNew.yJumps = imag(dataNew.yJumps);
            
        end
        
        % Insert a NaN (or array of NaNs) and append new data to array:
        xNaN = NaN(1, size(dataNew.xLine, 2)); % Array of NaNs.
        yNaN = NaN(1, size(dataNew.yLine, 2)); % Array of NaNs.

        data.xLine = [data.xLine ; xNaN ; dataNew.xLine];
        data.yLine = [data.yLine ; yNaN ; dataNew.yLine];
        data.xPoints = [data.xPoints ; xNaN ; dataNew.xPoints];
        data.yPoints = [data.yPoints ; yNaN ; dataNew.yPoints];
        data.xJumps = [data.xJumps ; dataNew.xJumps];
        data.yJumps = [data.yJumps ; dataNew.yJumps];

        data.xDeltas = [data.xDeltas ; dataNew.xDeltas];
        data.yDeltas = [data.yDeltas ; dataNew.yDeltas];
        data.yDeltaBase = [data.yDeltaBase ; dataNew.yDeltaBase];
        
        data.defaultXLim = data.defaultXLim & dataNew.defaultXLim;
        data.defaultYLim = data.defaultYLim & dataNew.defaultYLim;
        
    end

elseif ( nargin == 2 )
    % PLOT(F, G)
    
    [f, g] = overlap(f, g);

    % Loop over each FUN for Line and Points data:
    nFuns = numel(f.funs);
    
    for k = 1:nFuns
        % Get the data from the FUN objects:
        dataNew = plotData(f.funs{k}, g.funs{k});

        if ( ~any(ismember(fields(dataNew), 'xDeltas' )) )
            dataNew.xDeltas = [];
            dataNew.yDeltas = [];
            dataNew.yDeltaBase = [];
        end

        data.xLim = [min(dataNew.xLim(1), data.xLim(1)), ...
            max(dataNew.xLim(2), data.xLim(2))];
        data.yLim = [min(dataNew.yLim(1), data.yLim(1)), ...
            max(dataNew.yLim(2), data.yLim(2))];
        
        % Discard the unnecessary jump data:
        if ( k == 1 )
            dataNew.xJumps(1,:) = [];
            dataNew.yJumps(1,:) = [];
        end
        
        if ( k == nFuns )
            dataNew.xJumps(end,:) = [];
            dataNew.yJumps(end,:) = [];
        end

        % Array of NaNs:
        xNaN = NaN(1, size(dataNew.xLine, 2)); 
        yNaN = NaN(1, size(dataNew.yLine, 2)); 
        
        % Insert a NaN (or array of NaNs) and append new data to array:
        data.xLine = [data.xLine ; xNaN ; dataNew.xLine];
        data.yLine = [data.yLine ; yNaN ; dataNew.yLine];
        data.xPoints = [data.xPoints ; xNaN ; dataNew.xPoints];
        data.yPoints = [data.yPoints ; yNaN ; dataNew.yPoints];
        data.xJumps = [data.xJumps ; dataNew.xJumps];
        data.yJumps = [data.yJumps ; dataNew.yJumps];
        data.xDeltas = [data.xDeltas; dataNew.xDeltas];
        data.yDeltas = [data.yDeltas; dataNew.yDeltas];
        data.yDeltaBase = [data.yDeltaBase; dataNew.yDeltaBase];
        
        data.defaultXLim = data.defaultXLim & dataNew.defaultXLim;
        data.defaultYLim = data.defaultYLim & dataNew.defaultYLim;
    end
    
else
    % PLOT(F, G, H)
    
    % Initialise z storage:
    data.zLine = [];
    data.zPoints = [];
    data.zJumps = [];

    % Overlap the CHEBFUN objects:
    [f, g] = overlap(f, g);
    [g, h] = overlap(g, h);
    [h, f] = overlap(h, f);
    
    % Loop over each FUN for Line and Points data:
    nFuns = numel(f.funs);
    for k = 1:nFuns
        % Get the data from the FUN objects:
        dataNew = plotData(f.funs{k}, g.funs{k}, h.funs{k});

        if ( ~any(ismember(fields(dataNew), 'xDeltas' )) )
            dataNew.xDeltas = [];
            dataNew.yDeltas = [];
        end

        data.xLim = [min(dataNew.xLim(1), data.xLim(1)), ...
            max(dataNew.xLim(2), data.xLim(2))];
        data.yLim = [min(dataNew.yLim(1), data.yLim(1)), ...
            max(dataNew.yLim(2), data.yLim(2))];

        myNaN = NaN(1, size(dataNew.yLine, 2)); % Array of NaNs.
        % Insert a NaN (or array of NaNs) and append new data to array:
        data.xLine = [data.xLine ; myNaN ; dataNew.xLine];
        data.yLine = [data.yLine ; myNaN ; dataNew.yLine];
        data.zLine = [data.zLine ; myNaN ; dataNew.zLine];
        data.xPoints = [data.xPoints ; myNaN ; dataNew.xPoints];
        data.yPoints = [data.yPoints ; myNaN ; dataNew.yPoints];
        data.zPoints = [data.zPoints ; myNaN ; dataNew.zPoints];
        
    end
    
    % Return NaNs if there are no jumps:
    data.xJumps = myNaN;
    data.yJumps = myNaN;
    data.zJumps = myNaN;
    
    % Loop over each FUN for Jumps data:
    for k = 1:(nFuns - 1)
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

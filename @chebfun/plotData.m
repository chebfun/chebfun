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
    'xJumps', [], 'yJumps', [], 'yLim', []);

if ( nargin == 1 )
    % PLOT(F)
    
    % Overhead:
    nFuns = numel(f.funs);
    ymax = zeros(1,nFuns);
    ymin = zeros(1,nFuns);

    % Loop over each FUN for Line and Points data:
    for k = 1:nFuns
        % Get the data from the FUN:
        dataNew = plotData(f.funs{k});
        
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
        
        % If any of the boundary values is positively infinite, then set the
        % ylim for this piece to be the minimum of 10 and the largest yLine 
        % data:
        ymax(k) = max(max(dataNew.yLine));
        
        if ( dataNew.yLim(2) == Inf )
            ymax(k) = min(10, ymax(k));
        end
        
        % If any of the boundary values is negatively infinite, then set the
        % ylim for this piece to be the maximum of -10 and the smallest yLine 
        % data:
        ymin(k) = min(min(dataNew.yLine));
        
        if ( dataNew.yLim(1) == -Inf )
            ymin(k) = max(-10, ymin(k)); 
        end

    end
    
    % Take the maximum ymax to be the upper ylim for the entire CHEBFUN, while 
    % take the minimum ymin to be the lower ylim. Then store yLim in data:
    data.yLim = [min(ymin) max(ymax)];

elseif ( nargin == 2 )
    % PLOT(F, G)
    
    [f, g] = overlap(f, g);

    % Overhead:
    nFuns = numel(f.funs);
    ymax = zeros(1,nFuns);
    ymin = zeros(1,nFuns);
    
    for k = 1:nFuns
        % Get the data from the FUN objects:
        dataNew = plotData(f.funs{k}, g.funs{k});
        
        % Discard the unnecessary jump data:
        if ( k == 1 )
            dataNew.xJumps(1) = [];
            dataNew.yJumps(1,:) = [];
        end
        
        if ( k == nFuns )
            dataNew.xJumps(end) = [];
            dataNew.yJumps(end,:) = [];
        end

        % Array of NaNs:
        myNaN = NaN(1, size(dataNew.yLine, 2)); 
        
        % Insert a NaN (or array of NaNs) and append new data to array:
        data.xLine = [data.xLine ; myNaN ; dataNew.xLine];
        data.yLine = [data.yLine ; myNaN ; dataNew.yLine];
        data.xPoints = [data.xPoints ; myNaN ; dataNew.xPoints];
        data.yPoints = [data.yPoints ; myNaN ; dataNew.yPoints];
        data.xJumps = [data.xJumps ; dataNew.xJumps];
        data.yJumps = [data.yJumps ; dataNew.yJumps];
    end
    
    % If any of the boundary values is positively infinite, then set the
    % ylim for this piece to be 10:
    ymax(k) = max(max(dataNew.yLine));
    
    if ( dataNew.yLim(2) == Inf )
        ymax(k) = min(10, ymax(k));
    end
    
    % If any of the boundary values is negatively infinite, then set the
    % ylim for this piece to be -10:
    ymin(k) = min(min(dataNew.yLine));
    
    if ( dataNew.yLim(1) == -Inf )
        ymin(k) = max(-10, ymin(k));
    end
    
    % Store ylim:
    data.yLim = [min(ymin) max(ymax)];
    
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
    nFuns = numel(f.funs);
    for k = 1:nFuns
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
    
    % Return NaNs if there are no jumps:
    data.xJumps = NaN;
    data.yJumps = myNaN;
    
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

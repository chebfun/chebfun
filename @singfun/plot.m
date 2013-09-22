function varargout = plot(f, varargin)
%PLOT   Basic linear plot for SINGFUN objects. 
%   PLOT(F) plots the SINGFUN object F.
%
% See also PLOTDATA.
%
% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Deal with an empty input:
if ( isempty(f) )
    if ( nargout == 1 )
        varargout{1} = plot([]);
    end
    return
end

%%
% Store the hold state of the current axis:
holdState = ishold;

%%
% Get the data from SINGFUN.PLOTDATA():
data = plotData(f);

%%
% Plot the curve:
if ( isreal(data.fLine) )
    h1 = plot(data.xLine, data.fLine, varargin{:}); 
else
    h1 = plot(data.fLine, varargin{:}); 
end
set(h1, 'Marker', 'none') 
hold on

%%
% Plot the points:
if ( isreal(data.fLine) )
    h2 = plot(data.xPoints, data.fPoints, varargin{:});
else
    h2 = plot(data.fPoints, varargin{:});
end

%%
% Change the style accordingly:
set(h2,'LineStyle', 'none')
if ( all(strcmp(get(h2, 'Marker'), 'none')) && (length(f) < 257) )
    set(h2, 'Marker', 'o')
end

%% 
% [TODO]: This cell may be deleted. Added here for convenience.
% Since a SINGFUN usually blows up, set the y-axis limits to [-10,10]?
ylim( [-10, 10] )

%%
% Return hold state to what it was before:
if ( ~holdState )
    hold off
end

%%
% Give an output if one was requested:
if ( nargout > 0 )
    varargout{1} = h1;
    varargout{2} = h2;
end

end
function varargout = area(f,varargin)
% AREA   Filled CHEBFUN area plot.
%
% See also PLOT.

% Copyright 2013 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% Deal with empty case:
if ( isempty(f) )
    varargout{1} = plot(); %#ok<LTARG>
    return
end

% Get the (x, y) data from PLOTDATA():
data = plotData(f);
x = data.xLine; 
y = data.fLine;

% Remove NaNs from jumps:
if ( isnan(x(1)) )
    x(1) = x(2);
    y(1,:) = y(2,:);
end
if ( isnan(x(end)) )
    x(end) = x(end-1);
    y(end,:) = y(end-1,:);
end
mask = find(isnan(x));
% Take averae of left/right for y:
y(mask,:) = (y(mask+1,:) + y(mask-1,:))/2;
% Take subsequent point for x:
x(mask) = x(mask+1);

% Call built-in AREA():
h = area(x, y, varargin{:});

% Output handle to plot window:
if ( nargout == 1 )
    varargout{1} = h;
end

end

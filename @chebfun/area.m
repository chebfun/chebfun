function varargout = area(f, varargin)
%AREA   Filled CHEBFUN area plot.
%   AREA(X, F) or AREA(F) is the same as PLOT(X, F) or PLOT(F) except that the
%   area between 0 and F is filled. When F is a array-valued, AREA(F) plots the
%   columns of Y as filled areas.
%
%   AREA(F, LEVEL) specifies the base level for the area plot to be at the value
%   y = LEVEL. The default value is LEVEL = 0.
%
%   AREA(..., 'Prop1', VALUE1, 'Prop2', VALUE2,...) sets the specified
%   properties of the underlying areaseries objects.
% 
%   AREA(AX, ...) plots into axes with handle AX. Use GCA to get the handle to
%   the current axes or to create one if none exist.
% 
%   H = AREA(...) returns a vector of handles to areaseries objects.
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
if ( (nargin > 1) && isa(varargin{1}, 'chebfun') )
    g = varargin{1};
    varargin(1) = [];
    data = plotData(f, g);
else
    data = plotData(f);
end
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
% Take average of left/right for y:
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

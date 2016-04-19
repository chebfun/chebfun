function varargout = waterfall(varargin)
%WATERFALL   Waterfall plot of a SPHEREFUN.
%   WATERFALL(F) displays the waterfall plot of F.
%
%   WATERFALL(F, S) displays the column and row chebfuns of F that are used for
%   its approximation.  This is a 3D version of plot(f,S), where S is a string
%   (see PLOT).
%
%   WATERFALL(F, S, 'nslices', N) displays the first min(N,length(f)) column
%   and rows.
%
%   WATERFALL supports passing options to the plot, similar to standard Matlab
%   plot commands. The options supported are:
%       'color':      Color of lines and markers plotted.
%       'marker':     Marker for pivot points.
%       'markersize': Size of markers plotted at pivot points.
%
%   H = WATERFALL(...) returns a handle to a waterfall plot object.
%
% See also SPHEREFUN/PLOT.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

[varargout{1:nargout}] = waterfall@separableApprox(varargin{:});

end

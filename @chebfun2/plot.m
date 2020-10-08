function varargout = plot(varargin)
%PLOT  Surface plot of a CHEBFUN2.
%
%   PLOT(F) gives a surface plot of the CHEBFUN2 F, the same as SURF(F).
%   If F is complex-valued, it gives a phase portrait.
%
%   PLOT(F, 'zebra') gives a "zebra plot", black for values < 0
%   and white for values >= 0.
%
%   PLOT(F, S) Plotting with option string plots the column and row slices, and
%   pivot locations used in the construction of F.
%
%   When the first argument in options is a string giving details about
%   linestyle, markerstyle or colour, pivot locations are plotted. Various
%   line types, plot symbols and colors may be obtained with plot(F,S) where S
%   is a character string made from one element from any or all the following 3
%   columns, similar as in the usual plot command:
%
%           b     blue          .     point              -     solid
%           g     green         o     circle             :     dotted
%           r     red           x     x-mark             --    dashed
%           c     cyan          +     plus               -.    dashdot
%           m     magenta       *     star             (none)  no line
%           y     yellow        s     square
%           k     black         d     diamond
%                               v     triangle (down)
%                               ^     triangle (up)
%                               <     triangle (left)
%                               >     triangle (right)
%                               p     pentagram
%                               h     hexagram
%
%   For phase portraits see: E. Wegert, Visual Complex Functions: An
%   Introduction with Phase Portraits, Springer Basel, 2012, or for MATLAB code
%   to produce many different styles of phase portraits go to
%   http://www.visual.wegert.com
%
% See also SURF, MESH, PHASEPLOT.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

[varargout{1:nargout}] = plot@separableApprox(varargin{:});

end

function varargout = plot(varargin)
%PLOT  Surface plot of a CHEBFUN2.
%
%   PLOT(F) if F is a real-valued CHEBFUN2 then this is the surface plot and is
%   the same as surf(F). If F is a complex valued then this returns a domain
%   colouring plot of F.
%
%   PLOT(F) if F is a complex-valued CHEBFUN2 then we do Wegert's phase portrait
%   plots.
%
%   PLOT(F, S) Plotting with option string plots the column and row slices, and
%   pivot locations used in the construction of F.
%
%   When the first argument in options is a string giving details about
%   linestyle, markerstyle or colour then pivot locations are plotted. Various
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
%   introduction with Phase Portraits, Springer Basel, 2012, or for MATLAB code
%   to produce many different styles of phase portraits go to:
%   http://www.visual.wegert.com
%
% See also SURF, MESH.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

[varargout{1:nargout}] = plot@separableApprox(varargin{:});

end

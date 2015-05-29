function varargout = plot( f, varargin )
%PLOT  Surface plot of a SPHEREFUN.
%
%   PLOT(F) if F is a real-valued SPHEREFUN then this is the surface plot and is
%   the same as surf(F). If F is a complex valued then this returns a domain
%   colouring plot of F.
%
%   PLOT(F) if F is a complex-valued SPHEREFUN then we do Wegert's phase portrait
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

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information. 

if ( ~isempty(varargin) )
    domOld = f.domain;
    % Have the domain include the doubled sphere
    f.domain = f.domain - [0 0 pi 0];
end

h = plot@separableApprox( f, varargin{:} );

if ( ~isempty(varargin) )
    holdState = ishold;
    if ~holdState
        hold on;
    end
    % Plot the "doubled" sphere box.
        
    % Get domain and change to include doubled sphere
    dom = domOld - [0 0 pi pi];
    x = dom(1); 
    y = dom(3);
    w = dom(2)-dom(1);
    ht = dom(4)-dom(3);
    LW = 'LineWidth';
    LS = 'LineStyle';

    % Calculate boundary of plotting region.
    sw = dom(1) - w/4;
    se = dom(2) +  w/4;
    nw = dom(3) - ht/4;
    ne = dom(4) +  ht/4 + pi;
    
    % Remove the current rectangle and plot our own.
    delete(findobj(gca,'Type','Rectangle'));
    
    rectangle('Position', [x y w ht], LW, 2, LS, ':');
    rectangle('Position', [x y+pi w ht], LW, 2, LS, '-');
    axis( [sw se nw ne] ); axis equal;
    if ~holdState
        hold off;
    end
end

if ( nargout > 0 )
    varargout = { h }; 
end

end

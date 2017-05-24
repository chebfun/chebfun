function varargout = plot( f, varargin )
%PLOT  Surface plot of a SEPARABLEAPPROX.
%
%   PLOT(F) gives a surface plot of the SEPARABLEAPPROX F, the same
%   as SURF(F).  If F is complex-valued, it gives a phase portrait.
%
%   PLOT(F, 'zebra') gives a "zebra plot", black for values < 0
%   and white for values >= 0.
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

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information. 

holdState = ishold;
if ( ~isempty(varargin) )
    % See if first set of option makes it a pivot plot.
    
    if ( length(varargin{1}) < 5 )
    
        %% Column, row, pivot plot

        % Only option with <=3 letters is a colour, marker, line
        ll = regexp( varargin{1}, '[-:.]+', 'match' );
        cc = regexp( varargin{1}, '[bgrcmykw]', 'match' );       % color
        mm = regexp( varargin{1}, '[.ox+*sdv^<>ph]', 'match' );  % marker
        
        if ( ~isempty(ll) )
            if ( strcmpi(ll{1},'.') )
                % so we have . first. Don't plot a line.
                ll = {};
            elseif ( strcmpi(ll{1},'.-') )
                % so we have . first. Don't plot a line.
                ll{1} = '-';
            end
        end
        plotline = ~isempty(ll);  % plot row and col pivot lines?
        if ( isempty(mm) )
            mm{1}= '.';
        end
        if ( isempty(cc) ) 
            cc{1}= 'b';
        end
        % Plot the crosses and pivots on domain.
        
        % Get domain.
        dom = f.domain;
        x = dom(1); y = dom(3);
        w = dom(2)-dom(1);
        ht = dom(4)-dom(3);
        LW = 'LineWidth';
        
        % Calculate boundary of plotting region.
        sw = dom(1) - w/4;
        se = dom(2) +  w/4;
        nw = dom(3) - ht/4;
        ne = dom(4) +  ht/4;
        
        % Plot the square domain:
        rectangle('Position', [x y w ht], LW, 2); hold on;
        axis( [sw se nw ne] ); axis equal;
        
        % Plot pivots:
        crosses = f.pivotLocations;
        defaultopts = { 'MarkerSize', 7 };
        extraopts = { 'Marker', mm{:}, 'LineStyle', 'none', 'Color', cc{:} };
        if ( length(varargin) < 2 )
            h = plot( crosses(:,1), crosses(:,2), extraopts{:}, defaultopts{:} );
        else
            h = plot(crosses(:,1), crosses(:,2), extraopts{:},...
                                      defaultopts{:}, varargin{2:end} );
        end
        if ( plotline )
            opts = {}; 
            if ( ~isempty(ll) )
                opts = { opts{:}, 'LineStyle', ll{:}, 'color', cc{:} };
            end
            % Plot column lines:
            line( [crosses(:,1) crosses(:,1)].', [nw+ht/8 ne-ht/8], opts{:} );
            % Plot row lines:
            line( [sw+w/8 se-w/8], [crosses(:,2), crosses(:,2)].', opts{:} );
        end
        hold off
    else
        
        %% Standard or zebra surface plot 
        if strcmp(varargin{1}, 'zebra')
            h = surf(f);
            view(0,90)
            caxis(norm(caxis,inf)*[-1 1])
            colormap([0 0 0; 1 1 1])
        else
            h = surf(f, varargin{:});
        end
    end
else
    if ( isreal( f ) )
        h = surf( f );
        colormap default
        
    else
        %% Phase Protrait plot 
        % The following is a slightly modified version of Wegert's code from
        % p345 of his book.
        dom = f.domain; 
        nx = 500; 
        ny = 500;
        x = linspace(dom(1), dom(2), nx); 
        y = linspace(dom(3), dom(4), ny);
        [xx, yy] = meshgrid(x, y); 
        zz = xx + 1i*yy;
        f = feval(f, xx, yy);
        h = surf( real(zz), imag(zz), 0*zz, angle(-f) );
        set(h, 'EdgeColor', 'none');
        caxis([-pi pi]), colormap hsv(600)
        view(0, 90), axis equal, axis off  
        
    end
end

if ( ~holdState )
    hold off
end
if ( nargout > 0 )
    varargout = { h }; 
end

end

function varargout = plot( f, varargin )
%PLOT  Surface plot of a DISKFUN.
%
%   PLOT(F) gives a surface plot of the DISKFUN F, the same as SURF(F).
%   If F is complex-valued, it gives a phase portrait.
%
%   PLOT(F, 'zebra') gives a "zebra plot", black for values < 0
%   and white for values >= 0.
%
%   PLOT(F, S) with option string plots the column and row slices, and 
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
%
% See also DISKFUN/SURF.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Make a user option?
plot_full_grid = false;

if ( ~isempty(varargin) )
    
    % Plot the pivots on the surface of the disk.
    if ( length(varargin{1}) < 5 )
        dom = f.domain;
        holdState = ishold;
        if ( ~holdState )
            hold on
        end
        N = 100; % Used to generate solid disk as well as plot slicing lines
        th = trigpts(N, dom); th=[th; dom(2)];
        % If the plot is not being added to another then plot a unit circle
        % so the disk is more easily discernable.
        if ( ~holdState )
            % Generate a unit disk
            N = 200;
            th = trigpts(N, dom(1:2)); th=[th; dom(2)];
            r = exp(1i*th);
            plot(real(r),imag(r),'--', 'Linewidth', .5); 
        end
                
        % 
        % Column, row, pivot plot
        %
        
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
        plotline = ~isempty(ll);  % Plot row and col pivot lines?
        if ( isempty(mm) )
            mm{1}= '.';
        end
        if ( isempty(cc) )
            cc{1}= 'b';
        end
        % Plot the crosses and pivots on domain.
        
        % Full grid on the disk
        m = length(f.cols);
        n = length(f.rows);
        r = chebpts(m);
        r = r(floor(m/2)+1:m);
        [TT, RR] = meshgrid([trigpts(max(n,1000),dom(1:2)); dom(2)],r);
        [TT1, RR1] = meshgrid([trigpts(n,dom(1:2)); dom(2)],r);
        
        
        % Plot pivots:
        % Convert pivots to Cartesian coordinates
        pivots = f.pivotLocations;
        pivotsCart = zeros(size(pivots,1),2);
        pivotsCart(:,1) = (pivots(:,2)).*cos(pivots(:,1)); %x=rcosth
        pivotsCart(:,2) = (pivots(:,2)).*sin(pivots(:,1)); %y=rsinth
        XX = RR.*cos(TT);
        YY = RR.*sin(TT);
        XX1 = RR1.*cos(TT1);
        YY1 = RR1.*sin(TT1);

        if ( plot_full_grid )
            % Plot grayed out grid
            clrGrid = [192 192 192]/256;
            plot(XX1,YY1,'-','Color',clrGrid);
            plot(XX',YY','-','Color',clrGrid);
        end            
      
        % Also plot points marking the pivots shifted by pi in longitude
        % since these are technically also included in the GE algorithm.
        pivotsCart = [pivotsCart;-pivotsCart(:,1:2)];
        
        defaultopts = { 'MarkerSize', 7 };
        extraopts = { 'Marker', mm{:}, 'LineStyle', 'none', 'Color', cc{:} };
        if ( length(varargin) < 2 )
            h = plot(pivotsCart(:,1), pivotsCart(:,2), extraopts{:}, ...
                defaultopts{:});
        else
            h = plot(pivotsCart(:,1), pivotsCart(:,2), extraopts{:}, ...
                defaultopts{:}, varargin{2:end});
        end
        if ( plotline )

            colSlices = [];
            rowCircs = [];
            
            for k=1:size(pivots,1)
                if ( abs(pivotsCart(k,1)) > 100*eps )  %special case if x=0
                    temp = [cos(pivots(k,1)) sin(pivots(k,1))];
                    colSlices = [colSlices; temp;-temp; [nan nan]];
                    %rowcircs use theta pts set up earlier
                    if ( pivotsCart(k,2) == 0 )
                        rowCircs = [rowCircs; 0];
                        rowCircs = [rowCircs;nan];
                    else
                        rowCircs = [rowCircs; pivots(k,2)*exp(1i*th)];
                        rowCircs = [rowCircs; nan];
                    end
                else %case where x=0 so slope is undef
                    colSlices = [colSlices; zeros(2,1) [-1 ; 1]];
                    colSlices = [colSlices; nan nan];
                    if ( pivotsCart(k, 2) == 0 )
                        rowCircs = [rowCircs; 0];
                        rowCircs = [rowCircs; nan];
                    else
                        rowCircs = [rowCircs; pivots(k,2)*exp(1i*th)];
                        rowCircs = [rowCircs; nan];
                    end
                end
            end
            opts = {};
            if ( ~isempty(ll) )
                opts = { opts{:}, 'LineStyle', ll{:}, 'color', cc{:} };
            end
            
            % Plot column lines:
            plot(colSlices(:,1), colSlices(:,2), opts{:});
            
            % Plot row lines:
            plot(rowCircs, opts{:});
        end
        if ( ~holdState )
            hold off
        end
    else

        %% Standard or zebra surface plot 
        if strcmp(varargin{1}, 'zebra')
            h = surf(f);
            caxis(norm(caxis,inf)*[-1 1])
            colormap([0 0 0; 1 1 1])
       
            % Add unit circle to zebra plot
            holdState = ishold; 
            circ = exp(1i*pi*linspace(-1,1,101));
            hold on 
            plot(real(circ), imag(circ), 'k-','Color',0.2*[1 1 1])
            if ~holdState
                hold off;
            end
        else
            h = surf(f, varargin{:});
        end
    end
else
    h = plot@separableApprox(f);
end

if ( nargout > 0 )
    varargout = { h };
end
axis square
%axis off
view(2)
end

function varargout = plot(f, varargin)
%PLOT  Surface plot of a SPHEREFUN.
%
%   PLOT(F) if F is a real-valued SPHEREFUN then this is the surface plot 
%   and is the same as surf(F). If F is complex valued then this returns a 
%   domain colouring plot of F.
%
%   PLOT(F) if F is a complex-valued SPHEREFUN then we do Wegert's phase 
%   portrait plots.
%
%   PLOT(F, S) Plotting with option string plots the column and row slices,
%   and pivot locations used in the construction of F.
%
%   When the first argument in options is a string giving details about
%   linestyle, markerstyle or colour then pivot locations are plotted. 
%   Various line types, plot symbols and colors may be obtained with 
%   plot(F,S) where S is a character string made from one element from any 
%   or all the following 3 columns, similar as in the usual plot command:
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
%   Introduction with Phase Portraits, Springer Basel, 2012, or for MATLAB 
%   code to produce many different styles of phase portraits go to:
%   http://www.visual.wegert.com
% 
% See also SURF, MESH.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information. 

% Make a user option?
plot_full_grid = false;

if ( ~isempty(varargin) )
    
    % Plot the pivots on the surface of the sphere.
    if ( length(varargin{1}) < 5 )
        dom = f.domain;
        holdState = ishold;
        if ~holdState
            hold on
        end
        
        % If the plot is not being added to another then plot a solid 
        % sphere so the lines are more easily discernable.
        if ( ~holdState )
            % Generate a unit sphere.
            [XX,YY,ZZ] = sphere(101);
            
            % Color of the sphere will be yellow:
            clr = [255 255 204]/255;
            
            % Plot the sphere, make it slightly smaller than unit so lines
            % show up more clearly.
            scl = 0.99;
            surf(scl*XX, scl*YY, scl*ZZ, 1+0*XX, 'EdgeColor', 'None',...
                'FaceColor', clr);
        end

        %% Column, row, pivot plot
        
        % Only option with <=3 letters is a colour, marker, line
        ll = regexp(varargin{1}, '[-:.]+', 'match');
        cc = regexp(varargin{1}, '[bgrcmykw]', 'match');       % color
        mm = regexp(varargin{1}, '[.ox+*sdv^<>ph]', 'match');  % marker
        
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

        % Full grid on the sphere
        [n, m] = length(f);
        [LL,TT] = meshgrid(linspace(dom(1), dom(2), n+1), ...
            linspace(dom(3), dom(4), m/2+1));

        % Plot pivots:
        
        % Convert pivots to Cartesian coordinates
        pivots = f.pivotLocations;
        pivotsCart = zeros(size(pivots, 1), 3);
        if ( iscolat(f) )
            pivotsCart(:, 1) = cos(pivots(:, 1)).*sin(pivots(:, 2));
            pivotsCart(:, 2) = sin(pivots(:, 1)).*sin(pivots(:, 2));
            pivotsCart(:, 3) = cos(pivots(:, 2));
            XX = cos(LL).*sin(TT);
            YY = sin(LL).*sin(TT);
            ZZ = cos(TT);
        else
            pivotsCart(:, 1) = cos(pivots(:, 1)).*cos(pivots(:, 2));
            pivotsCart(:, 2) = sin(pivots(:, 1)).*cos(pivots(:, 2));
            pivotsCart(:, 3) = sin(pivots(:, 2));
            XX = cos(LL).*cos(TT);
            YY = sin(LL).*cos(TT);
            ZZ = sin(TT);
        end
        if ( plot_full_grid )
            % Plot grayed out grid
            clrGrid = [192 192 192]/256;
            plot3(XX, YY, ZZ, '-', 'Color', clrGrid);
            plot3(XX', YY', ZZ', '-', 'Color', clrGrid);
        end
        % Also plot points marking the pivots shifted by pi in longitude
        % since these are technically also included in the GE algorithm.
        pivotsCart = [ pivotsCart; [ -pivotsCart(:,1:2) pivotsCart(:,3) ] ];
            
        defaultopts = {'MarkerSize', 7};
        extraopts = { 'Marker', mm{:}, 'LineStyle', 'none', 'Color', cc{:} };
        if ( length(varargin) < 2 )
            h = plot3( pivotsCart(:,1), pivotsCart(:,2), ...
                       pivotsCart(:,3), extraopts{:}, defaultopts{:} );
        else
            h = plot3(pivotsCart(:, 1), pivotsCart(:, 2), ...
                      pivotsCart(:, 3), extraopts{:}, defaultopts{:}, ...
                      varargin{2:end});
        end
        if ( plotline )
            % Use parametrization for great circles passing through each
            % pivot location and the poles for the column pivot lines and
            % horizontal circles at the right latitude for the row pivot
            % lines.
            if ( plot_full_grid )
                num_t = max(50*size(XX, 1), 201);
            else
                num_t = 201;
            end
                
            t = linspace(-pi, pi, num_t)';  % Parametric variable
            colCircs = [];
            rowCircs = [];
            baseRowCirc = [ cos(t) sin(t) ];
            for k=1:size(pivots, 1)
                % Need to deal with the poles separately
                if ( (abs(pivotsCart(k, 3) + 1) > 100*eps) && ...
                        (abs(pivotsCart(k, 3) -1) > 100*eps) )
                    r = pivotsCart(k, 1:2)/sqrt(sum(pivotsCart(k, 1:2).^2));
                    colCircs = [ colCircs; sin(t)*r cos(t) ];
                    colCircs = [ colCircs; nan nan nan ];
                    rowCircs = [ rowCircs; sin(pivots(k, 2))*baseRowCirc ...
                        pivotsCart(k, 3)*ones(num_t, 1)];
                    rowCircs = [ rowCircs; nan nan nan];
                else % Pivot at the pole
                    if ( ~isempty(colCircs) ) % No circle through poles yet
                        colCircs = [ colCircs; sin(t)*[1 0 0] cos(t) ];
                        colCircs = [ colCircs; nan nan nan ];
                    end
                end
            end
            opts = {}; 
            if ( ~isempty(ll) )
                opts = { opts{:}, 'LineStyle', ll{:}, 'color', cc{:} };
            end
            % Plot column lines:
            plot3(colCircs(:, 1), colCircs(:, 2), colCircs(:, 3), ...
                opts{:}, 'linewidth', 1);
            % Plot row lines:
            plot3(rowCircs(:, 1), rowCircs(:, 2), rowCircs(:, 3), ...
                opts{:},'linewidth', 1);            
        end

        if ( ~holdState )
            hold off
            view(180-37.5, 30);
            daspect([1 1 1]);
        end
    else
        %% Standard surface plot 
        h = surf(f, varargin{:});
    end
else
    h = plot@separableApprox(f);
end

if ( nargout > 0 )
    varargout = { h }; 
end

end
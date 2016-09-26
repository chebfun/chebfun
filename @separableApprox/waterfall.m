function varargout = waterfall( f, varargin )
%WATERFALL   Waterfall plot of a SEPARABLEAPPROX.
%   WATERFALL(F) displays the waterfall plot of F.
%
%   WATERFALL(F, S) displays the column and row chebfuns of F that are used for
%   its approximation.  This is a 3D version of plot(f,S), where S is a string
%   (see PLOT).
%
%   WATERFALL(F, S, 'nslices', N) displays the first min(N,length(f)) columns
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
% See also PLOT.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

numpts = 200;

% Empty check: 
if ( isempty(f) ) 
    h = plot( [] );
    if ( nargin > 0 )
        varargout = { h };
    end
end

% Options for plotting:
plotOpts = {};

% For plotting, it's useful to know whether we're running in old or new
% Matlab graphics mode
if ( ~verLessThan('matlab', '8.4') )
    newMatlabVersion = true;
else
    newMatlabVersion = false;
end

% Number of points to slices:
nslices = length(f); 
j = 1; argin = {};
while ( ~isempty(varargin) )
    if ( strcmpi(varargin{1}, 'nslices') )
        nslices = varargin{2};
        varargin(1:2) = [];
    elseif ( strcmpi(varargin{1}, 'color') )
        plotOpts = [plotOpts, 'color', varargin{2}];
        varargin(1:2) = [];
    elseif ( strcmpi(varargin{1}, 'marker') )
        plotOpts = [plotOpts, 'marker', varargin{2}];
        varargin(1:2) = [];
    elseif ( strcmpi(varargin{1}, 'markersize') )
        plotOpts = [plotOpts, 'markersize', varargin{2}];
        varargin(1:2) = [];
    else
        argin{j} = varargin{1};
        varargin(1) = [];
        j = j+1;
    end
end

nslices = min( nslices, length(f) );
varargin = argin; 
mynargin = length(varargin) + 1; 

% Get domain:
dom = f.domain;     

if ( mynargin == 1 )     
    % Standard waterfall plot
    
    x = linspace(dom(1), dom(2), numpts);
    y = linspace(dom(3), dom(4), numpts);
    [xx, yy] = meshgrid(x, y);
    vals = feval( f, xx, yy );
    h = waterfall( xx, yy, vals );
    return
    
else
    % See if first set of option makes it a pivot plot.
    
    if ( length(varargin{1})<5 )
        % Column and row waterfall plot
        
        holdState = ishold;

        % Only option with <= 3 letters is a colour, marker, line
        ll = regexp(varargin{1},'[-:.]+','match');
        cc = regexp(varargin{1},'[bgrcmykw]','match');        % color
        mm = regexp(varargin{1},'[.ox+*sdv^<>ph]', 'match');  % marker
        
        if ( ~isempty(ll) )
            if ( strcmpi(ll{1},'.') )
                % so we have . first. Don't plot a line.
                ll = {};
            elseif ( strcmpi(ll{1},'.-') )
                % so we have . first. Don't plot a line.
                ll{1} = '-';
            end
        end
        % Plot row and col pivot lines?
        plotline = ~isempty(ll); 
        if ( isempty(mm) ) 
            mm{1}= '.';
        end
        if ( isempty(cc) )
            cc{1}= 'b';
        end
        % Evaluate on a grid.
        P = f.pivotLocations;
        
        defaultopts = {'markersize', 20, 'marker', '.'};
        lineopts = {'linewidth', 2, 'linestyle', ll{:}};
        if ( ~plotline )
            % Just plot the pivots at height f(x,y)
            h = plot3(P(:,1), P(:,2), feval( f, P(:,1), P(:,2) ), '.', ...
                defaultopts{:}, plotOpts{:});
        else
            h = [];
%             hold on
%             % Plot the pivots, rows and columns, one set at a time:
%             for pivotCounter = 1:nslices
%                 % Reset color cycle prior to point plot if running R2014b.
%                 if ( newMatlabVersion )
%                     set(gca, 'ColorOrderIndex', pivotCounter);
%                 end
%                 tempPoint = P(pivotCounter,:)
%                 hTemp = plot3(tempPoint(1), tempPoint(2), ...
%                     feval(f, tempPoint(1), tempPoint(2)), ...
%                     defaultopts{:}, plotOpts{:});
% %                 h1 = plot3(P(pivotCounter, 1), P(ss, 2), feval( f, P(ss ,1), P(ss, 2) ), ...
% %                     '.', defaultopts{:}, plotOpts{:} ); hold on
%                 pivotCounter
%                 h = [h; hTemp];
%             end

            ss = 1:nslices;
            
            
            [xx, yy]=meshgrid (P(:, 1), chebpts( length(f.cols), dom(3:4) ) );
            vals = feval( f, xx, yy );
            
            % Plot column slices
            for jj = 1:nslices
                xxTemp = chebfun( P(jj, 1), dom(1:2) );
                yyTemp = chebfun( [-1 ; 1], dom(3:4) ); 
                zzTemp = chebfun( vals(:,jj), dom(1:2) );
                if ( newMatlabVersion )
                    set(gca, 'ColorOrderIndex', jj);
                end
                plot3(P(jj, 1), P(jj, 2), feval( f, P(jj ,1), P(jj, 2) ), ...
                    '.', defaultopts{:}, plotOpts{:} ); hold on
                if ( newMatlabVersion )
                    set(gca, 'ColorOrderIndex', jj);
                end
                h2 = plot3( xxTemp, yyTemp, zzTemp, lineopts{:} ); hold on
            end
            
            [xx, yy] = meshgrid( chebpts( length(f.rows), dom(1:2) ), P(:,2) );
            vals = feval( f, xx, yy );            
            % Plot row slices:
            for jj = 1:nslices
                xxTemp = chebfun( [-1 ; 1], dom(1:2) ); 
                yyTemp = chebfun( P(jj, 2), dom(3:4) );
                zzTemp = chebfun( vals(jj,:).', dom(3:4) );
                if ( newMatlabVersion )
                    set(gca, 'ColorOrderIndex', jj);
                end
                h3 = plot3( xxTemp, yyTemp, zzTemp, lineopts{:} );
            end
            axis equal
            
%             h = [h1(:) ; h2(:) ; h3(:)];
            if ( ~holdState )
                hold off
            end
        end
    end
end

%% Output

if ( nargout > 0 )
    varargout = {h};
end

end

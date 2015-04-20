function varargout = waterfall( f, varargin )
%WATERFALL   Waterfall plot of a CHEBFUN2.
%   WATERFALL(F) displays the waterfall plot of F.
%
%   WATERFALL(F, S) displays the column and row chebfuns of F that are used for
%   its approximation.  This is a 3D version of plot(f,S), where S is a string
%   (see PLOT).
%
%   WATERFALL(F, S, 'nslices', N) displays the first min(N,length(f)) columns
%   and rows.
%
%   H = WATERFALL(...) returns a handle to a waterfall plot object.
%
% See also PLOT.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

numpts = 200;

% Empty check: 
if ( isempty(f) ) 
    h = plot( [] );
    if ( nargin > 0 )
        varargout = { h };
    end
end

% Number of points to slices:
nslices = length(f); 
j = 1; argin = {};
while ( ~isempty(varargin) )
    if strcmpi(varargin{1}, 'nslices')
        nslices = varargin{2};
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
        
        defaultopts = {'MarkerSize', 7};
        extraopts = {'Marker', mm{:}, 'Color', cc{:}};
        lineopts = {'linewidth', 2, 'linestyle', ll{:}};
        if ( ~plotline )
            % Just plot the pivots at height f(x,y)
            h = plot3(P(:,1), P(:,2), feval( f, P(:,1), P(:,2) ), '.', ...
                                             extraopts{:}, defaultopts{:});
        else
            % Plot the pivots.
            ss = 1:nslices;
            h1 = plot3(P(ss, 1), P(ss, 2), feval( f, P(ss ,1), P(ss, 2) ),...
                                '.', extraopts{:}, defaultopts{:} ); hold on
            
            
            [xx, yy]=meshgrid (P(:, 1), chebpts( length(f.cols), dom(3:4) ) );
            vals = feval( f, xx, yy );
            
            % Plot column slices
            xx = [];
            yy = []; 
            zz = [];
            for jj = 1:nslices
                xx = [xx chebfun( P(jj, 1), dom(1:2) )];
                yy = [yy chebfun( [-1 ; 1], dom(3:4) )]; 
                zz = [zz chebfun( vals(:,jj), dom(1:2) )];
            end
            h2 = plot3( xx, yy, zz, lineopts{:} ); hold on
            
            [xx, yy] = meshgrid( chebpts( length(f.rows), dom(1:2) ), P(:,2) );
            vals = feval( f, xx, yy );
            xx = []; 
            yy = []; 
            zz = [];
            
            % Plot row slices:
            for jj = 1:nslices
                xx = [xx chebfun( [-1 ; 1], dom(1:2) )]; 
                yy = [yy chebfun( P(jj, 2), dom(3:4) )];
                zz = [zz chebfun( vals(jj,:).', dom(3:4) )];
            end
            h3 = plot3( xx, yy, zz, lineopts{:} );
            axis equal
            
            h = [h1(:) ; h2(:) ; h3(:)];
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

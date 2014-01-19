function varargout = surfc(f,varargin)
%SURFC  Combination surf/contour plot for a chebfun2.
% 
% SURFC(...) is the same as SURF(...) except that a contour plot
% is drawn beneath the surface.
% 
% See SURFC, SURFL.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.

if ( isempty( f ) ) % check for empty chebfun2.
    h=surfc([]);  % call the empty surfc command.
    if nargout == 1
        varargout = {h};
    end
    return;
end

minplotnum = 200; % How dense to make the samples.
defaultopts = {'facecolor','interp','edgealpha',.5,'edgecolor','none'};

% Number of points to plot
j = 1; argin = {};
while ( ~isempty(varargin) )
    if ( strcmpi(varargin{1}, 'numpts') )
        minplotnum = varargin{2};
        varargin( 1:2 ) = [];
    else
        argin{j} = varargin{1};
        varargin(1) = [];
        j = j+1;
    end
end

if ( isempty(argin) )
    argin = {};
end

if ( isa(f,'chebfun2') )                 % surfc(f,...)
    if ( ( nargin == 1 ) ||...
            ( nargin > 1 && ~isempty(argin) && ~isa(argin{1},'chebfun2') ) ||...
            ( nargin == 3 && isempty(argin))) 
        % Get domain.
        dom = f.domain;
        x = chebfun2(@(x,y) x, dom); y = chebfun2(@(x,y) y, dom);
        h = surfc(x, y, f,defaultopts{:}, argin{:}, 'numpts', minplotnum);
        xlim( dom(1:2) ), ylim( dom(3:4) )
    elseif ( nargin > 2)             %surfc(x,y,f,...), with x, y, f chebfun2 objects
        % Extract arguments: 
        x = f; 
        y = argin{1};
        if ( isa(y,'chebfun2') )
            % check domains of x and y are the same.
            if ( ~chebfun2.domainCheck(x, y) )
                error('CHEBFUN2:surfc:DATADOMAINS','Domains of chebfun2 objects do not match.');
            end
        end
        dom = x.domain; 
        xdata = linspace(dom(1), dom(2), minplotnum);
        ydata = linspace(dom(3), dom(4), minplotnum);
        [xx, yy] = meshgrid(xdata, ydata);
        x = feval(x, xx, yy); 
        y = feval(y, xx, yy);
        if ( isa(argin{2}, 'chebfun2') )      % surfc(x,y,f,...)
            vals = feval(argin{2}, xx, yy);
            if ( nargin < 4 )  % surfc(x,y,f)
                C = vals;
            elseif ( isa(argin{3},'double') )    % surfc(x,y,f,C,...)
                C = argin{3};
                argin(3)=[];
            elseif ( isa(argin{3},'chebfun2'))  % colour matrix given as a chebfun2.
                C = feval(argin{3},xx,yy);
                argin(3)=[];
            else
                C = vals;
            end
            
            % make some correct to C, for prettier plotting.
            if ( norm(C - C(1,1), inf) < 1e-10 )
                % If vals are very close up to round off then the color scale is
                % hugely distorted.  This fixes that.
                [n, m] = size( C );
                C = C(1, 1) * ones(n, m);
            end
            
            h = surfc( x, y, vals, C, defaultopts{:}, argin{3:end} );
            xlabel('x'), ylabel('y'), xlim(dom(1:2)), ylim(dom(3:4))
            
            % There is a bug in matlab surfc plot when vals are very nearly a constant.
            % Fix this manually by resetting axis scaling.
            if ( norm( vals - vals(1,1), inf) < 1e-10 * norm(vals, inf) )...
                     && ( ~(norm(vals - vals(1,1),inf)==0) )
                v = vals(1,1); absv = abs(v);
                zlim([v-.5*absv v+.5*absv])
            end
            
        else
            error('CHEBFUN2:surfc:INPUTS','The third argument should be a chebfun2 if you want to supply chebfun2 data.')
        end
    else                                            %surfc(f,C)
        dom = f.corners;
        x = chebfun2(@(x,y) x,dom); y = chebfun2(@(x,y) y,dom);
        h = surfc(x,y,f,argin{1},defaultopts{:},argin{2:end});
        xlim(dom(1:2)), ylim(dom(3:4))
    end
else     % surfc(X,Y,f,...)
    error('CHEBFUN2:surfc:INPUTS','Data should be given as chebfun2 objects \n For example: \n x = chebfun2(@(x,y)x); y = chebfun2(@(x,y)y);\n surfc(x,y,f)');
end


if ( nargout > 0 )
    varargout = {h};
end

end
function varargout = quiver3( F, varargin )
%QUIVER3   3-D quiver plot of a CHEBFUN2V.
%   QUIVER3(F) plots velocity vectors as arrows with components F(1), F(2),
%   F(3), which are CHEBFUN2 objects. QUIVER3 automatically scales the arrows to
%   fit. The arrows are plotted on a uniform grid.
%
%   QUIVER3(Z,F) plots velocity vectors at the equally spaced surface points
%   specified by the matrix or CHEBFUN2 Z. If Z is a CHEBFUN2 then we use Z to
%   map the uniform grid.
%
%   QUIVER3(X,Y,Z,F) plots velocity vectors at (x,y,z). If X, Y, Z are CHEBFUN2
%   objects then we use X, Y, Z to map the uniform grid.
%
%   QUIVER3(F,S), QUIVER3(Z,F,S) or QUIVER3(X,Y,Z,F,S) automatically scales the
%   arrows to fit and then stretches them by S. Use S=0 to plot the arrows with
%   the automatic scaling.
%
%   QUIVER3(...,LINESPEC) uses the plot linestyle specified for the velocity
%   vectors.  Any marker in LINESPEC is drawn at the base instead of an arrow on
%   the tip.  Use a marker of '.' to specify no marker at all.  See PLOT for
%   other possibilities.
%
%   QUIVER(...,'numpts',N) plots arrows on a N by N uniform grid.
%
%   QUIVER3(...,'filled') fills any markers specified.
%
%   H = QUIVER3(...) returns a quiver object.
%
%   If F is a CHEBFUN2V with two components then we recommend using
%   CHEBFUN2V/QUIVER.
%
% See also QUIVER.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

numpts = 20; 

% Empty check:
if ( isempty( F ) )
    quiver([])
    return
end

if ( isempty(varargin) )
    varargin = {};
end

% Number of points to plot
j = 1; 
argin = {};
while ( ~isempty(varargin) )
    if strcmpi(varargin{1}, 'numpts')
        numpts = varargin{2};
        varargin(1:2) = [];
    else
        argin{j} = varargin{1};
        varargin(1) = [];
        j = j+1;
    end
end
varargin = argin; 

if ( isa(F, 'chebfun2v') && ...
        (isempty(varargin) || ~isa(varargin{1}, 'chebfun2v')) ) 
    % quiver3(F,...)
    
    nF = F.nComponents; 
    if ( nF == 2 ) 
        h = quiver( F, varargin{:} ); 
    else
        % Plot arrows are equally spaced locations
        dom = F.components{1}.domain;
        x = linspace(dom(1), dom(2), numpts);
        y = linspace(dom(3), dom(4), numpts);
        [xx, yy] = meshgrid(x, y); 
        zz = zeros( size( xx ) );
        vals1 = feval(F.components{1}, xx, yy);
        vals2 = feval(F.components{2}, xx, yy);
        vals3 = feval(F.components{3}, xx, yy);
        h = quiver3( xx, yy, zz, vals1, vals2, vals3, varargin{:} );
    end

elseif ( isa(F, 'chebfun2') && isa(varargin{1}, 'chebfun2v') )  
    % quiver(Z,F,...)
    
    % Extract arguments:
    Z = F; 
    F = varargin{1};
    % Domain check:
    if ( ~domainCheck(Z, F.components{1} ) )
        error('CHEBFUN:CHEBFUN2V:quiver3:domain', ...
            'Object are not on the same domain.');
    end
    % Workout plotting locations:
    zz = new_data_locations( Z, numpts );
    %Plot: 
    h = quiver3( zz, F, varargin{2:end} );
    
elseif ( isa(F, 'double') && isa(varargin{1}, 'chebfun2v') )  
    % quiver(zz,F,...)
    
    nF = F.nComponents; 
    if ( nF == 2 ) 
        h = quiver( F, varargin{:} ); 
    else
        % Plot arrows are equally spaced locations
        dom = F.components{1}.domain;
        x = linspace(dom(1), dom(2), numpts);
        y = linspace(dom(3), dom(4), numpts);
        [xx, yy] = meshgrid(x, y); 
        vals1 = feval(F.components{1}, xx, yy);
        vals2 = feval(F.components{2}, xx, yy);
        vals3 = feval(F.components{3}, xx, yy);
        h = quiver3( F, vals1, vals2, vals3, varargin{:} );
    end
    
elseif ( nargin > 3 )        
    % quiver(xx,yy,zz,F,...) or quiver(X,Y,Z,F,...)
    
    if ( isa(F,'double') )       
        % quiver(xx,yy,zz,F,...)
        
        % Extract arguments: 
        xx = F; 
        yy = varargin{1}; 
        zz = varargin{2}; 
        F = varargin{3};
        % Check that we have the right input types: 
        if ( ~isa(yy,'double') || ~isa(zz,'double') || ~isa(F,'chebfun2v') )
            error('CHEBFUN:CHEBFUN2V:quiver3:badInputs1',...
                'Unrecognised input arguments.');
        end
        % Get plotting data: 
        dom = F.components{1}.domain;
        xdata = linspace( dom(1), dom(2), size(xx,1) );
        ydata = linspace (dom(3), dom(4), size(yy,2) ); 
        [mxdata, mydata]=meshgrid(xdata, ydata); 
        vals1 = feval(F.components{1}, mxdata, mydata);
        vals2 = feval(F.components{2}, mxdata, mydata);
        vals3 = feval(F.components{3}, mxdata, mydata);
        % Plot:
        h = quiver3(xx, yy, zz, vals1, vals2, vals3, varargin{4:end});
        
    elseif isa(F,'chebfun2')                         
        % quiver(X,Y,Z,F,...)
        
        % Extract arguments: 
        X = F; 
        Y = varargin{1}; 
        Z = varargin{2}; 
        F = varargin{3};
        % Check that we have the right input types: 
        if ( ~isa(Y,'chebfun2') || ~isa(Z,'chebfun2') || ~isa(F,'chebfun2v') )
            error('CHEBFUN:CHEBFUN2V:quiver3:badInputs2', ...
                'Unrecognised input arguments.');
        end
        % Check domains: 
        if ( ~domainCheck(X, Y) || (~domainCheck(Y, Z) ) )
            error('CHEBFUN:CHEBFUN2V:quiver3:domain',...
                'Object are not on the same domain.');
        end
        
        % Get new data locations.
        xx = new_data_locations(X, numpts);
        yy = new_data_locations(Y, numpts);
        zz = new_data_locations(Z, numpts);
        
        % Plot quiver3.
        h = quiver3( xx, yy, zz, F, varargin{4:end} );
    end
    
else
    error('CHEBFUN:CHEBFUN2V:quiver3:badInputs3',...
        'Unrecognised input arguments.');
end

if ( nargout > 0 )
    varargout = { h };
end

end

function newloc = new_data_locations( f1 , numpts )
% Generate new arrow location if first two inputs are CHEBFUN2 objects.

% Check the CHEBFUN2 objects are on the same domain.
dom = f1.domain;

% mesh 'em up for the quiver arrows.
x = linspace(dom(1), dom(2), numpts);
y = linspace(dom(3), dom(4), numpts);

[xx, yy] = meshgrid(x, y);
% Use CHEBFUN2 to generate data locations.
newloc = feval(f1, xx, yy);      

end

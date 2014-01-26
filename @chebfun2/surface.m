function varargout = surface( f, varargin )
%SURFACE  Plot surface of a chebfun2.
% 
% surface(X,Y,Z,C) adds the surface in X,Y,Z,C to the current axes.
% surface(X,Y,Z) uses C = Z, so color is proportional to surface height.
% See SURF for a complete description of the various forms that X,Y,Z,C
% can take.
% 
% See SURF. 

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.

if ( isempty( f ) )
     surface([]); 
     return
end

ish = ishold;
numpts = 200; 
dom = f.domain; 

% Data points: 
[xx, yy] = chebfun2.chebpts2( numpts, numpts, dom );
val = feval(f, xx, yy);

% options: 
defaultopts = {'facecolor','interp','edgealpha',.5,'edgecolor','none'};


if ( isempty(varargin) )
    h1 = surface(xx, yy, val, defaultopts{:} ); 
    hold on 
    h2 = surface( xx.', yy.', val.', defaultopts{:} );
else
    h1 = surface( xx, yy, val, defaultopts{:}, varargin{:} ); 
    hold on 
    h2 = surface(xx.', yy.', val.', defaultopts{:}, varargin{:} );
end

if ( ~ish ) 
    hold off 
end 

if ( nargout > 0 ) 
    varargout = {h1 h2}; 
    return
end

end
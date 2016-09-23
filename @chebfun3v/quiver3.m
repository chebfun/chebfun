function varargout = quiver3(F, varargin)
%QUIVER3   3-D quiver plot of a CHEBFUN3V object.
%   QUIVER3(F) plots velocity vectors as arrows with components F(1), F(2)
%   and F(3) which are CHEBFUN3 objects. QUIVER3 automatically scales the 
%   arrows to fit. The arrows are plotted on a uniform grid.
%
%   QUIVER3(F, S) automatically scales the arrows to fit and then stretches 
%   them by S. Use S=0 to plot the arrows with the automatic scaling.
%
%   QUIVER3(..., LINESPEC) uses the plot linestyle specified for the 
%   velocity vectors.  Any marker in LINESPEC is drawn at the base instead 
%   of an arrow on the tip. Use a marker of '.' to specify no marker at 
%   all.
%
% See also quiver3.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

numpts = 7; 

% Empty check:
if ( isempty(F) )
    h = quiver3([], [], [], []);
    return
end

if ( isa(F, 'chebfun3v') && ...
        (isempty(varargin) || ~isa(varargin{1}, 'chebfun3v')) ) 
    
    nF = F.nComponents;
    if ( nF == 2 )
        error('CHEBFUN:CHEBFUN3V:quiver3:badInputs1',...
            'Needs three CHEBFUN3 objects.');
    else
        % Plot arrows are equally spaced locations
        dom = F.components{1}.domain;
        x = linspace(dom(1), dom(2), numpts);
        y = linspace(dom(3), dom(4), numpts);
        z = linspace(dom(5), dom(6), numpts);
        [xx, yy, zz] = ndgrid(x, y, z);
        vals1 = feval(F.components{1}, xx, yy, zz);
        vals2 = feval(F.components{2}, xx, yy, zz);
        vals3 = feval(F.components{3}, xx, yy, zz);
        h = quiver3(xx, yy, zz, vals1, vals2, vals3, varargin{:});
        axis([x(1), x(end), y(1), y(end), z(1), z(end)])        
    end

else
    error('CHEBFUN:CHEBFUN3V:quiver3:badInputs3',...
        'Unrecognised input arguments.');
end

if ( nargout > 0 )
    varargout = {h};
end

end
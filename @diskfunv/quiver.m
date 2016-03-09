function varargout = quiver( F, varargin )
%QUIVER   Quiver plot of DISKFUNV.
%   QUIVER(F) plots the vector velocity field of F. QUIVER automatically
%   attempts to scale the arrows to fit within the grid. The arrows are on a
%   uniform grid.
%
% BELOW HERE NOT YET IMPLEMENTED: 
%   QUIVER(F,S) automatically scales the arrows to fit within the grid and then
%   stretches them by S.  Use S=0 to plot the arrows without the automatic
%   scaling. The arrows are on a uniform grid.
%
%   QUIVER(X,Y,F,...) is the same as QUIVER(F,...) except the arrows are on the
%   grid given in X and Y.
%
%   QUIVER(...,LINESPEC) uses the plot linestyle specified for the velocity
%   vectors.  Any marker in LINESPEC is drawn at the base instead of an arrow on
%   the tip.  Use a marker of '.' to specify no marker at all.  See PLOT for
%   other possibilities.
%
%   QUIVER(...,'numpts',N) plots arrows on a N by N uniform grid.
%
%   H = QUIVER(...) returns a quivergroup handle.
%
%   If F is a DISKFUN with three non-zero components then this calls
%   QUIVER3. (this is currently not implemented in diskfun)
%
% See also QUIVER3.
% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

numpts = 30;

% Empty check:
if ( isempty( F ) )
    quiver([])
    return
end

if ( isempty(varargin) )
    varargin = {};
end

dom = [-pi pi 0 1]; 

%add solid disk when there is no other plot in order to make arrows more visible.
holdState = ishold;
if ~holdState
    hold on;
end

        if ~holdState
            %
            % Generate a unit disk
            N = 200;
            th = trigpts(N, dom(1:2)); th=[th; dom(2)];
            r = exp(1i*th);
            
            clr = [255 255 204]/255;
            fill(real(r),imag(r), clr, 'Edgecolor', 'None');             
        end

% Number of points to plot
j = 1;
argin = {};
while ( ~isempty( varargin ) )
    if strcmpi( varargin{1}, 'numpts' )
        numpts = varargin{2};
        varargin(1:2) = [];
    else
        argin{j} = varargin{1};
        varargin(1) = [];
        j = j+1;
    end
end
varargin = argin;

if ( isa(F, 'diskfunv') )             % quiver(F,...)
    
    nF = F.nComponents;
    if ( nF == 3 )
        h = quiver3(F, varargin{:});   % Call quiver3 instead.
    else
        % Plot quiver with arrows at equally spaced points: 
        [xx, yy] = diskpts(numpts, 0, 0);
        F1 = F.components{1}; F2 = F.components{2};
        vals1 = feval(F1, xx, yy, 'cart'); 
        vals2 = feval(F2, xx, yy, 'cart');
        h = quiver(xx, yy, vals1,vals2, varargin{:});
        axis(1.1*[-1 1 -1 1])
    end
    
elseif ( nargin >= 3 )                 % quiver(x,y,F,...)
    
    % First two arguments contain arrow locations: %for diskfun we need to
    % know if this is cartesian or polar...for now we assume cartesian
    
    x = F;
    y = varargin{1};
    
    if ( isa(varargin{2}, 'diskfunv') )
        F = varargin{2};
        nF = F.nComponents;
        if ( nF == 3 )
            h = quiver3(F,varargin{:}); % Call quiver3 instead.
        else
            [xx, yy] = diskpts(0,x,y); 
            F1 = F.components{1}; F2 = F.components{2};
            vals1 = feval(F1, xx, yy, 'cart');
            vals2 = feval(F2, xx, yy, 'cart');
            %dom = F1.domain;
            h = quiver( xx, yy, vals1, vals2, varargin{3:end} );
            axis(1.1*[-1 1 -1 1]);
        end
    else
        error('DISKFUN:DISKFUNV:quiver:inputs', ...
            'Third argument should be a diskfunv.');
    end
    
end

if ( nargout > 0 )
    varargout = {h};
end


 function [xx,yy] = diskpts(s, x,y)  % gets good spaced points on disk in cartesian coords. Will need improvement.
                                    %for user-input coords, this throws out
                                    % choices that are not on the unit
                                    % disk. 
if (s > 0) 
x = linspace(-1,1,s); 
y = linspace(-1,1,s);
[xx, yy] = meshgrid(x,y); 
else
xx = x; 
yy = y;
end

xx=xx(:);
yy = yy(:);

R = xx.^2+yy.^2;   
a = find(bsxfun(@max, R, ones(length(xx),1))-R); %throw out stuff outside unit disk                                           
xx = xx(a);
yy = yy(a);

end
   

end

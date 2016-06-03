function varargout = scan(f, varargin)
%SCAN   Contour plots of a CHEBFUN3 object.
%   SCAN(F, DIM) plot contours of a CHEBFUN3 object F at 9 different
%   cross sections with the coordinate DIM, where DIM is either 1, 2 or 3.
%
%   SCAN(F, DIM, 'HOLD') plot contours of a CHEBFUN3 object F at 3 cross 
%   sections along the coordinate DIM and hold all three.
%
%   SCAN(F) is the same as SCAN(F, 1).
%
%   If F is complex-valued, then phase portraits are used.
% 
% See also CHEBFUN3/PLOT, CHEBFUN3/SLICE, CHEBFUN3/ISOSURFACE, and 
% CHEBFUN3/SURF.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

holdState = ishold;
noSlices = 9; 
noHoldSlices = 3;
pauseTime = 0.5;
dom = f.domain;
isHold = 0; % The default 'hold' flag.

if ( nargin == 1 )
    dim = 1; 
elseif ( nargin == 2 )
    if ( strcmpi(varargin{1}, 'hold'))
    isHold = 1; 
    dim = 1;
    elseif ( strcmpi(varargin{1}, 'ct'))
    ct_flag = 1; 
    dim = 1;    
    else
        dim = varargin{1};
    end
elseif ( nargin == 3 )
    if ( strcmpi(varargin{1}, 'hold') )
        isHold = 1; 
        dim = varargin{2};
    elseif ( strcmpi(varargin{2}, 'hold') )
        isHold = 1; 
        dim = varargin{1};
    end
end

%% Determine some values of the function
numpts = 51;
[xx,yy,zz] = meshgrid(linspace(dom(1), dom(2), numpts), ...
    linspace(dom(3), dom(4), numpts), linspace(dom(5), dom(6), numpts));
vv = feval(f, xx, yy, zz);

%%
if ( ( dim == 3 ) && ~isHold )
    zslices = linspace(dom(5), dom(6), noSlices);
    
    for i=1:numel(zslices)
        if isreal(vv)
            h = slice(xx, yy, zz, vv, [], [], zslices(i));
            set(h, 'EdgeColor', 'none');
            colorbar
        else
            h = slice(xx, yy, zz, angle(-vv), [], [], zslices(i));
            set(h, 'EdgeColor', 'none');
            caxis([-pi pi])
            colormap('hsv')
            axis('equal')
        end
        axis(dom)
        hold off
        title(['Slice: z = ' num2str(zslices(i))])
        pause(pauseTime)
        drawnow
    end
elseif ( ( dim == 3 ) && isHold )
    zslices = linspace(dom(5), dom(6), noHoldSlices);
    for i=1:numel(zslices)
        if isreal(vv)
            h = slice(xx, yy, zz, vv, [], [], zslices(i));
            colorbar
        else
            h = slice(xx, yy, zz, angle(-vv), [], [], zslices(i));
            caxis([-pi pi])
            colormap('hsv')
            axis('equal')            
        end
        axis(dom)
        set(h,'EdgeColor','none');
        hold on
        title(['Slice: z = ' num2str(zslices(i))])
        pause(pauseTime)
        drawnow
    end
        
elseif ( ( dim == 2 ) && ~isHold )
    yslices = linspace(dom(3), dom(4), noSlices);
    for i=1:numel(yslices)
        if isreal(vv)
            h = slice(xx, yy, zz, vv, [], yslices(i), []);
            colorbar
        else
            h = slice(xx, yy, zz, angle(-vv), [], yslices(i), []);
            caxis([-pi pi])
            colormap('hsv')
            axis('equal')
        end
        set(h, 'EdgeColor', 'none')
        hold off
        axis(dom);
        title(['Slice: y = ' num2str(yslices(i))])
        pause(pauseTime)
        drawnow
    end
elseif ( ( dim == 2 ) && isHold )
    yslices = linspace(dom(3), dom(4), noHoldSlices);
    for i=1:numel(yslices)
        if ( isreal(vv) )
            h = slice(xx, yy, zz, vv, [], yslices(i), []);
            colorbar
        else
            h = slice(xx, yy, zz, angle(-vv), [], yslices(i), []);
            caxis([-pi pi])
            colormap('hsv')
            axis('equal')            
        end
        axis(dom)
        set(h, 'EdgeColor', 'none')
        hold on
        title(['Slice: y = ' num2str(yslices(i))])
        pause(pauseTime)
        drawnow
    end
elseif ( ( dim == 1 ) && ~isHold )
    xslices = linspace(dom(1), dom(2), noSlices);
    for i=1:numel(xslices)
        if ( isreal(vv) )
            h = slice(xx, yy, zz, vv, xslices(i), [], []);
            colorbar
        else
            h = slice(xx, yy, zz, angle(-vv), xslices(i), [], []);
            caxis([-pi pi])
            colormap('hsv')
            axis('equal')            
        end
        set(h, 'EdgeColor', 'none');
        hold off
        axis(dom);
        title(['Slice: x = ' num2str(xslices(i))])
        pause(pauseTime)
        drawnow
    end
elseif ( ( dim == 1 ) && isHold )
    xslices = linspace(dom(1), dom(2), noHoldSlices);
    for i=1:numel(xslices)
        if ( isreal(vv) )
            h = slice(xx, yy, zz, vv, xslices(i), [], []);
            colorbar
        else
            h = slice(xx, yy, zz, angle(-vv), xslices(i), [], []);
            caxis([-pi pi])
            colormap('hsv')
            axis('equal')            
        end
        set(h, 'EdgeColor', 'none')
        axis(dom);
        title(['Slice: x = ' num2str(xslices(i))])
        pause(pauseTime)
        drawnow
        hold on
    end
end

if ( ~holdState )
    hold off
end
if ( nargout > 0 )
    varargout = {h};
end

end % End of function
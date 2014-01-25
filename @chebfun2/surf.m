function varargout = surf( f, varargin )
%SURF  Surface plot of a chebfun2.
%
% SURF(F, C) plots the colored parametric surface defined by F and the
% matrix C. The matrix C, defines the colouring of the surface.
%
% SURF(F) uses colors proportional to surface height.
%
% SURF(X,Y,F,...) is the same as SURF(F,...) when X and Y are chebfun2
% objects except X and Y supplies the plotting locations are  mapped by
% X and Y.
%
% SURF(...,'PropertyName',PropertyValue,...) sets the value of the
% specified surface property.  Multiple property values can be set
% with a single statement.
%
% SURF returns a handle to a surface plot object.
%
% See also PLOT, SURFC.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.

if ( isempty( f ) ) % check for empty chebfun2.
    h=surf([]);  % call the empty surf command.
    if nargout == 1
        varargout = {h};
    end
    return;
end

minplotnum = 200; % How dense to make the samples.
defaultopts = {'facecolor','interp','edgealpha',.5,'edgecolor','none'};

% Number of points to plot
j = 1; argin = {};
while ( ~isempty(varargin) )                % SURF(F, 'NUMPTS', N)
    if ( strcmpi(varargin{1}, 'numpts') )
        minplotnum = varargin{2};
        varargin( 1:2 ) = [];
        dom = f.domain;
        x = linspace(dom(1), dom(2), minplotnum);
        y = linspace(dom(3), dom(4), minplotnum);
        [xx, yy] = meshgrid(x, y); 
        vf = feval(f, xx, yy);
        h = surf(xx, yy, vf, defaultopts{ : }, varargin{ : }); 
        if ( nargout > 0 ) 
            varargout = { h };
        end
        return
    else
        argin{j} = varargin{1};
        varargin(1) = [];
        j = j+1;
    end
end

if ( isempty(argin) )
    argin = {};
end

% The goal of the if-else statements are to convert the user defined SURF to one
% that can be dealt with by MATLAB's SURF command. 

% SURF( XX, YY, F , ...) 
if ( isa(f, 'double') ) 
    if ( nargin >= 3 ) 
        xx = f; 
        yy = argin{1}; 
        f = argin{2}; 
        argin(1:2) = []; 
        if ( isa( f, 'chebfun2' ) )
            if ( all( size(xx) == size(yy) ) )
                vf = feval(f, xx, yy);
                h = surf(xx, yy, vf, defaultopts{:}, argin{:});
            else
                error('CHEBFUN2:SURF:SIZES', 'Evaluation points size inconsistency.')
            end
        else
           error('CHEBFUN2:SURF:DISCRETE','Third argument should be a chebfun2.') 
        end
    else
        error('CHEBFUN2:SURF:INPUTS', 'Unrecognized input syntax.');
    end
elseif ( isa(f, 'chebfun2') )                 % SURF(F, ...)
    dom = f.domain;
        x = linspace(dom(1), dom(2), minplotnum);
        y = linspace(dom(3), dom(4), minplotnum);
        [xx, yy] = meshgrid(x, y); 
        vf = feval(f, xx, yy);
    if ( nargin == 1 )                        % SURF(F)      
        h = surf(xx, yy, vf, defaultopts{:}, argin{:}); 
    elseif ( nargin == 2)                     % SURF(F, C), C = chebfun2
        if ( isa(argin{1}, 'chebfun2') )
            vC = feval(argin{1}, xx, yy); 
            argin{1} = [];
            h = surf(xx, yy, vf, vC, defaultopts{:}); 
        else
            error('CHEBFUN2:SURF:INPUTS', 'Unrecognized input syntax.')
        end
    elseif ( nargin >= 3 )                    % SURF(X, Y, F)
        % Peel off user inputs:
        xcheb = f; 
        ycheb = argin{1}; 
        f = argin{2}; 
        if ( ~isa(xcheb, 'chebfun2') ||  ~isa(ycheb, 'chebfun2')...
                                                    ||  ~isa(f, 'chebfun2'))
            error('CHEBFUN2:SURF:INPUTS','Data should be given as chebfun2 objects.')
        end
        argin(1:2) = []; 
        if ( ~chebfun2.domainCheck(xcheb, ycheb) )
            error('CHEBFUN2:SURF:INPUTS','Inconsistent domains.');
        end
        if ( ~chebfun2.domainCheck(xcheb, f) )
            error('CHEBFUN2:SURF:INPUTS','Inconsistent domains.'); 
        end
        dom = xcheb.domain;
        x = linspace(dom(1), dom(2), minplotnum);
        y = linspace(dom(3), dom(4), minplotnum);
        [xx, yy] = meshgrid(x, y); 
        vf = feval(f, xx, yy);
        h = surf(xx, yy, vf, defaultopts{:}, argin{:}); 
    end
else
    error('CHEBFUN2:SURF:INPUTS', ...
        'Data should be given as chebfun2 objects, see help SURF');
end

if ( nargout > 0 )
    varargout = { h };
end

end
function varargout = surfc( f, varargin )
%SURFC  Combination surf/contour plot for a chebfun2.
% 
% SURFC(...) is the same as SURF(...) except that a contour plot
% is drawn beneath the surface.
% 
% See SURFC, SURFL.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.

if ( isempty( f ) ) % check for empty chebfun2.
    h=surfc([]);  % call the empty surf command.
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
        h = surfc(xx, yy, vf, defaultopts{ : }, varargin{ : }); 
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
                h = surfc(xx, yy, vf, defaultopts{:}, argin{:});
            else
                error('CHEBFUN2:SURFC:SIZES', 'Evaluation points size inconsistency.')
            end
        else
           error('CHEBFUN2:SURFC:DISCRETE','Third argument should be a chebfun2.') 
        end
    else
        error('CHEBFUN2:SURFC:INPUTS', 'Unrecognized input syntax.');
    end
elseif ( isa(f, 'chebfun2') )                 % SURF(F, ...)
    dom = f.domain;
        x = linspace(dom(1), dom(2), minplotnum);
        y = linspace(dom(3), dom(4), minplotnum);
        [xx, yy] = meshgrid(x, y); 
        vf = feval(f, xx, yy);
    if ( nargin == 1 )                        % SURF(F)      
        h = surfc(xx, yy, vf, defaultopts{:}, argin{:}); 
    elseif ( nargin == 2)                     % SURF(F, C), C = chebfun2
        if ( isa(argin{1}, 'chebfun2') )
            vC = feval(argin{1}, xx, yy); 
            argin{1} = [];
            h = surfc(xx, yy, vf, vC, defaultopts{:}); 
        else
            error('CHEBFUN2:SURFC:INPUTS', 'Unrecognized input syntax.')
        end
    elseif ( nargin >= 3 )                    % SURF(X, Y, F)
        % Peel off user inputs:
        xcheb = f; 
        ycheb = argin{1}; 
        f = argin{2}; 
        if ( ~isa(xcheb, 'chebfun2') ||  ~isa(ycheb, 'chebfun2')...
                                                    ||  ~isa(f, 'chebfun2'))
            error('CHEBFUN2:SURFC:INPUTS','Data should be given as chebfun2 objects.')
        end
        argin(1:2) = []; 
        if ( ~chebfun2.domainCheck(xcheb, ycheb) )
            error('CHEBFUN2:SURFC:INPUTS','Inconsistent domains.');
        end
        if ( ~chebfun2.domainCheck(xcheb, f) )
            error('CHEBFUN2:SURFC:INPUTS','Inconsistent domains.'); 
        end
        dom = xcheb.domain;
        x = linspace(dom(1), dom(2), minplotnum);
        y = linspace(dom(3), dom(4), minplotnum);
        [xx, yy] = meshgrid(x, y); 
        vf = feval(f, xx, yy);
        h = surfc(xx, yy, vf, defaultopts{:}, argin{:}); 
    end
else
    error('CHEBFUN2:SURFC:INPUTS', ...
        'Data should be given as chebfun2 objects, see help SURFC');
end

if ( nargout > 0 )
    varargout = { h };
end

end
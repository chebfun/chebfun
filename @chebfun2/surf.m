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

% if ( isempty( f ) ) % check for empty chebfun2.
%     h=surf([]);  % call the empty surf command.
%     if nargout == 1
%         varargout = {h};
%     end
%     return;
% end
% 
% minplotnum = 200; % How dense to make the samples.
% defaultopts = {'facecolor','interp','edgealpha',.5,'edgecolor','none'};
% 
% % Number of points to plot
% j = 1; argin = {};
% while ( ~isempty(varargin) )                % SURF(F, 'NUMPTS', N)
%     if ( strcmpi(varargin{1}, 'numpts') )
%         minplotnum = varargin{2};
%         varargin( 1:2 ) = [];
%         dom = f.domain;
%         x = linspace(dom(1), dom(2), minplotnum);
%         y = linspace(dom(3), dom(4), minplotnum);
%         [xx, yy] = meshgrid(x, y); 
%         vf = feval(f, xx, yy);
%         h = surf(xx, yy, vf, defaultopts{ : }, varargin{ : }); 
%         if ( nargout > 0 ) 
%             varargout = { h };
%         end
%         return
%     else
%         argin{j} = varargin{1};
%         varargin(1) = [];
%         j = j+1;
%     end
% end
% 
% if ( isempty(argin) )
%     argin = {};
% end
% 
% % The goal of the if-else statements are to convert the user defined SURF to one
% % that can be dealt with by MATLAB's SURF command. 
% 
% % SURF( XX, YY, F , ...) 
% if ( isa(f, 'double') ) 
%     if ( nargin >= 3 ) 
%         xx = f; 
%         yy = argin{1}; 
%         f = argin{2}; 
%         argin(1:2) = []; 
%         if ( isa( f, 'chebfun2' ) )
%             if ( all( size(xx) == size(yy) ) )
%                 vf = feval(f, xx, yy);
%                 h = surf(xx, yy, vf, defaultopts{:}, argin{:});
%             else
%                 error('CHEBFUN2:SURF:SIZES', 'Evaluation points size inconsistency.')
%             end
%         else
%            error('CHEBFUN2:SURF:DISCRETE','Third argument should be a chebfun2.') 
%         end
%     else
%         error('CHEBFUN2:SURF:INPUTS', 'Unrecognized input syntax.');
%     end
% elseif ( isa(f, 'chebfun2') )                 % SURF(F, ...)
%     dom = f.domain;
%         x = linspace(dom(1), dom(2), minplotnum);
%         y = linspace(dom(3), dom(4), minplotnum);
%         [xx, yy] = meshgrid(x, y); 
%         vf = feval(f, xx, yy);
%     if ( nargin == 1 )                        % SURF(F)      
%         h = surf(xx, yy, vf, defaultopts{:}, argin{:}); 
%     elseif ( nargin == 2)                     % SURF(F, C), C = chebfun2
%         if ( isa(argin{1}, 'chebfun2') )
%             vC = feval(argin{1}, xx, yy); 
%             argin{1} = [];
%             h = surf(xx, yy, vf, vC, defaultopts{:}); 
%         else
%             error('CHEBFUN2:SURF:INPUTS', 'Unrecognized input syntax.')
%         end
%     elseif ( nargin == 3 )                    % SURF(X, Y, F)
%         % Peel off user inputs:
%         xcheb = f; 
%         ycheb = argin{1}; 
%         f = argin{2}; 
%         if ( ~isa(xcheb, 'chebfun2') ||  ~isa(ycheb, 'chebfun2')...
%                                                     ||  ~isa(f, 'chebfun2'))
%             error('CHEBFUN2:SURF:INPUTS','Data should be given as chebfun2 objects.')
%         end
%         argin(1:2) = []; 
%         if ( ~chebfun2.domainCheck(xcheb, ycheb) )
%             error('CHEBFUN2:SURF:INPUTS','Inconsistent domains.');
%         end
%         if ( ~chebfun2.domainCheck(xcheb, f) )
%             error('CHEBFUN2:SURF:INPUTS','Inconsistent domains.'); 
%         end
%         dom = xcheb.domain;
%         x = linspace(dom(1), dom(2), minplotnum);
%         y = linspace(dom(3), dom(4), minplotnum);
%         [xx, yy] = meshgrid(x, y); 
%         vx = feval(xcheb, xx, yy); 
%         vy = feval(ycheb, xx, yy); 
%         h = surf(vx, vy, vf, defaultopts{:}, argin{:}); 
%     elseif ( nargin >= 4 )                    % SURF(X, Y, F)
%         % Peel off user inputs:
%         xcheb = f; 
%         ycheb = argin{1}; 
%         f = argin{2};
%         argin(1:2) = []; 
%         C = 0; 
%         if ( isa(argin{1}, 'chebfun2' ) )
%             C = argin{1};
%             argin(1) = []; 
%         end
%         if ( ~isa(xcheb, 'chebfun2') ||  ~isa(ycheb, 'chebfun2')...
%                                                     ||  ~isa(f, 'chebfun2'))
%             error('CHEBFUN2:SURF:INPUTS','Data should be given as chebfun2 objects.')
%         end
%         
%         if ( ~chebfun2.domainCheck(xcheb, ycheb) )
%             error('CHEBFUN2:SURF:INPUTS','Inconsistent domains.');
%         end
%         if ( ~chebfun2.domainCheck(xcheb, f) )
%             error('CHEBFUN2:SURF:INPUTS','Inconsistent domains.'); 
%         end
%         if ( isa(C,'chebfun2') &&  ~chebfun2.domainCheck(xcheb, f) )
%             error('CHEBFUN2:SURF:INPUTS','Inconsistent domains.'); 
%         end
%         dom = xcheb.domain;
%         x = linspace(dom(1), dom(2), minplotnum);
%         y = linspace(dom(3), dom(4), minplotnum);
%         [xx, yy] = meshgrid(x, y); 
%         vx = feval(xcheb, xx, yy); 
%         vy = feval(ycheb, xx, yy); 
%         if ( isa(C, 'chebfun2') )
%             vC = feval(C, xx, yy);
%             h = surf(vx, vy, vf, vC, defaultopts{:}, argin{:});
%         else
%             h = surf(vx, vy, vf, defaultopts{:}, argin{:});
%         end
%     end
% else
%     error('CHEBFUN2:SURF:INPUTS', ...
%         'Data should be given as chebfun2 objects, see help SURF');
% end
% 
% if ( nargout > 0 )
%     varargout = { h };
% end
% 
% end


if ( isempty(f) ) % check for empty chebfun2.
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
while ( ~isempty(varargin) )
    if strcmpi(varargin{1},'numpts')
        minplotnum = varargin{2};
        varargin(1:2) = [];
    else
        argin{j} = varargin{1};
        varargin(1) = [];
        j = j+1;
    end
end

if isempty(argin)
    argin = {};
end

if ( isa(f,'chebfun2') )
    if ( ( nargin == 1 ) || ( nargin > 1 && ~isempty(argin) && ~isa(argin{1},'chebfun2') ) || ( nargin == 3 && isempty(argin))) % surf(f,...)
        % Get domain.
        dom = f.domain;
        x = chebfun2(@(x,y) x,dom); y = chebfun2(@(x,y) y,dom);
        h = surf(x,y,f,defaultopts{:},argin{:},'numpts',minplotnum);
        xlim(dom(1:2)), ylim(dom(3:4))
    elseif ( nargin > 2)                    %surf(x,y,f,...), with x, y, f chebfun2 objects
        x = f; y = argin{1};
        if isa(y,'chebfun2')
            % check domains of x and y are the same.
            dom = x.domain; rectcheck = y.domain;
            if any(dom - rectcheck)
                error('CHEBFUN2:SURF:DATADOMAINS','Domains of chebfun2 objects do not match.');
            end
        end
        xdata = linspace(dom(1),dom(2),minplotnum);
        ydata = linspace(dom(3),dom(4),minplotnum);
        [xx,yy] = meshgrid(xdata,ydata);
        x = feval(x,xx,yy); 
        y = feval(y,xx,yy);
        if ( isa(argin{2},'chebfun2') )      % surf(x,y,f,...)
            vals = feval(argin{2},xx,yy);
            if nargin < 4   % surf(x,y,f)
                C = vals;
            elseif ( isa(argin{3},'double') )    % surf(x,y,f,C,...)
                C = argin{3};
                argin(3)=[];
            elseif ( isa(argin{3},'chebfun2'))  % colour matrix given as a chebfun2.
                C = feval(argin{3},xx,yy);
                argin(3)=[];
            else
                C = vals;
            end
            
            % make some correct to C, for prettier plotting.
            if ( norm(C - C(1,1),inf) < 1e-10 )
                % If vals are very close up to round off then the color scale is
                % hugely distorted.  This fixes that.
                [n,m]=size(C);
                C = C(1,1)*ones(n,m);
            end
            
            h = surf(x,y,vals,C,defaultopts{:},argin{3:end});
            xlabel('x'), ylabel('y'), xlim(dom(1:2)), ylim(dom(3:4))
            
            % There is a bug in matlab surf plot when vals are very nearly a constant.
            % Fix this manually by resetting axis scaling.
            if norm(vals - vals(1,1),inf)<1e-10*norm(vals,inf) && ~(norm(vals - vals(1,1),inf)==0)
                v = vals(1,1); absv = abs(v);
                zlim([v-.5*absv v+.5*absv])
            end
            
        else
            error('CHEBFUN2:SURF:INPUTS','The third argument should be a chebfun2 if you want to supply chebfun2 data.')
        end
    else  %surf(f,C)
        dom = f.domain;
        x = chebfun2(@(x,y) x,dom); y = chebfun2(@(x,y) y,dom);
        h = surf(x,y,f,argin{1},defaultopts{:},argin{2:end});
        xlim(dom(1:2)), ylim(dom(3:4))
    end
else     % surf(X,Y,f,...)
    error('CHEBFUN2:SURF:INPUTS','Data should be given as chebfun2 objects \n For example: \n x = chebfun2(@(x,y)x); y = chebfun2(@(x,y)y);\n surf(x,y,f)');
end


if nargout > 0
    varargout = {h};
end

end
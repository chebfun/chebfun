%DISKFUNV   Class constructor for DISKFUNV objects.
%
%   DISKFUNV(F, G) constructs a DISKFUNV with two components from the
%   function handles F and G. It represents vector-valued functions defined
%   on the disk of radius one. These can also be DISKFUN objects. Each
%   component of a DISKFUN represented as a DISKFUN.
%
% See also DISKFUN.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

classdef diskfunv
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CLASS PROPERTIES:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    properties ( Access = public )
        components      % Array of DISKFUN objects.
        isTransposed    % Has it been transposed?
        domain          % Domain of F
        nComponents     % Number of components.
        coords          % Is it in polar or Cartesian coordinates?
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CLASS METHODS:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods
        function f = set.coords(f, propName)
            if (strcmp(propName, 'polar') || strcmp(propName, 'cart'))
                f.coords = propName;
                f.components{1}.coords=propName;
                f.components{2}.coords=propName;
            else  %error if unacceptable setting is provided
                error('CHEBFUN:DISKFUNV:setcoords:propName', ...
                    'Coordinate setting must be either ''polar'' or ''cart''')
            end
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CLASS CONSTRUCTOR:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods ( Access = public, Static = false )
        
        function F = diskfunv( varargin )
            % DISKFUNV  The main DISKFUNV constructor!
            %
            % DISKFUNV is a class for representing vector-valued functions
            % on the disk. Each component is represented as a DISKFUN
            % object. This class allows one to conveniently compute with
            % vector-calculus on the disk.
            %
            % For more information see:
            %
            % A. Townsend, G. Wright, H. Wilber, "Computing with function
            % in polar and spherical geometries, II. The disk, submitted to
            % SISC, 2016.
            %
            % See also SPHEREFUNV
            
            % Return an empty DISKFUNV:
            if ( (nargin == 0) || isempty(varargin) )
                return
            end
            
            % This function calls the DISKFUN constructor once for each
            % non-zero component because a DISKFUNV is a vector of
            % DISKFUN objects.
            
            % If the argument is a DISKFUNV, nothing to do:
            if ( isa(varargin{1}, 'diskfunv') )
                F = varargin{1};
                return
            end
            
            % Stop now if there are too few input arguments:
            if ( nargin < 2 )
                error('CHEBFUN:DISKFUNV:diskfunv:arrayValued', ...
                    'Less than two components is not supported.')
            end
            
            %The domain will always be the following: (in polar coords)
            dom = [-pi pi 0 1];
            
            % Go find the vectorize flag in the user inputs (if there):
            vectorize = 0;
            for jj = 1:numel(varargin)
                if ( strcmpi( varargin{jj}, 'vectorize' ) )
                    vectorize = 1;
                    varargin(jj) = [];
                    break
                end
            end
            
            % Determine if the coordinate system has been provided by the
            % user and if it conflicts with coord systems. We default to
            % Cartesian coordinates.
            if isa(varargin{1}, 'diskfun') && isa(varargin{2}, 'diskfun')
                if ((nargin < 3) && (~strcmpi(varargin{1}.coords, varargin{2}.coords)) )
                    warning('CHEBFUN:DISKFUNV:diskfunv:coords',...
                        ['The two components have different coordinate',...
                        'settings. Now setting the diskfunv to evaluate',...
                        'with the default Cartesian coordinates.'])
                    iscart = 1;
                else
                    % If two diskfuns are supplied, then grap the coord
                    % setting:
                    iscart = diskfun.coordsetting(varargin{:});
                end
                
                % Diskfun and a function given; take diskfun setting and look
                % for coordinate flag
            elseif ( isa(varargin{1}, 'diskfun') &&...
                                     isa(varargin{2}, 'function_handle') )
                  
                % Coord system is set by the diskfun object: 
                iscart = diskfun.coordsetting(varargin);
                
            elseif ( isa(varargin{1}, 'function_handle') &&...
                                               isa(varargin{2}, 'diskfun') )
                
                % Coord system is set by the diskfun object:           
                iscart = diskfun.coordsetting(varargin{2}, varargin{3:end});
                
            else % Just go and check the coordinate flag. 
                
                ispolar = find(strcmp(varargin,'polar'));
                iscart = ~any(ispolar); 
                
            end
            
            % Go and set the coordset variable for the coordinary system to
            % construct the diskfunv object: 
            if ( ~iscart )
                coordset = 'polar';
            else
                coordset = 'cart';
            end
            
            % Get rid of coordinate argument in the user inputs so we 
            % can determine the number of components in our diskfunv object:
            for jj = 1:numel(varargin)
                if ( strcmpi( varargin{jj}, 'cart' ) ||...
                                         strcmpi( varargin{jj}, 'polar' ) )
                    varargin(jj) = [];
                    break
                end
            end
            
            % Unwrap input arguments; first two inputs are components
            component = 1;
            for jj = 1:numel(varargin)
                if ( iscell( varargin{jj} ) )
                    for kk = 1:numel( varargin{jj} )
                        fh{component} = varargin{jj}{kk};
                        component = component + 1;
                    end
                else
                    fh{component} = varargin{jj};
                    component = component + 1;
                end
            end
            varargin = fh;
            
            % DISKFUNV objects are vector-valued so complain if there is
            % only one component: 
            if ( numel(fh) < 2 )
                error('CHEBFUN:DISKFUNV:diskfunv:arrayValued', ...
                    'Less than two components is not supported.')
            end
            
            % DISKFUNV objects cannot contain more than two components. 
            % Complain if we have been given three or more.  
            if ( numel( fh ) > 2 )
                error('CHEBFUN:DISKFUNV:diskfunv:arrayValued', ...
                    'More than two components is not supported.')
            end
            
            % Now, convert all function handles for the DISKFUNV components
            % to DISKFUN objects:
            for jj = 1:2
                if ( isa( varargin{jj}, 'function_handle') )
                    % Function handle:
                    if ( ~vectorize )
                        newcheb = diskfun( varargin{jj}, coordset);
                    else
                        newcheb = diskfun( varargin{jj}, 'vectorize', coordset);
                    end
                    fh{jj} = newcheb;
                    
                elseif ( isa( varargin{jj}, 'diskfun') )
                    % DISKFUN object:
                    fh{jj} = varargin{jj};
                    fh{jj}.coords = coordset;
                    
                elseif ( isa( varargin{jj}, 'chebfun') )
                    % CHEBFUN object:
                    fh{jj} = diskfun( varargin{jj}, coordset);
                    
                elseif ( isa( varargin{jj}, 'double') )
                    % DOUBLE object:
                    fh{jj} = diskfun( varargin{jj}, coordset);
                    
                end
            end
            
            % Finally, assign the two components to make a DISKFUNV object:
            F.components = fh;
            F.isTransposed = 0;
            F.domain = dom;
            F.nComponents = numel( fh );
            F.coords = coordset;
            
        end
    end 
end
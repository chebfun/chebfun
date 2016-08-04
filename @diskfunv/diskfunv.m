%DISKFUNV   Class constructor for DISKFUNV objects.
% 
% DISKFUNV(F,G) constructs a DISKFUNV with two components from the
% function handles F and G.  These can also be DISKFUN objects or any
% other object that the DISKFUN constructor accepts.  Each component is
% represented as a DISKFUN.
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
        isTransposed    % transposed?
        domain          % Domain of F
        nComponents     % number of components. 
        coords          % polar or Cartesian coordinates.
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
    %% STATIC METHODS:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods ( Access = public, Static = true )
        % Unit normal vector
%         n = unormal( dom );   
        % Unit normal vector for r (radius) and theta (angle)
%         n = unit(a); 
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CLASS CONSTRUCTOR:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods ( Access = public, Static = false )
        
        function F = diskfunv( varargin )
            % The main DISKFUNV constructor!
                       
            % Return an empty DISKFUNV:
            if ( (nargin == 0) || isempty(varargin) )
                return
            end
                       
            % This function calls the DISKFUN constructor once for each 
            % non-zero component because a DISKFUNV  is just vector of 
            % DISKFUN objects.
            
            % If the argument is a DISKFUNV, nothing to do:
            if ( isa(varargin{1}, 'diskfunv') ) 
                F = varargin{1};
                return
            end
            
             % Stop now if there are too few components: 
            if ( nargin < 2 ) 
                error('CHEBFUN:DISKFUNV:diskfunv:arrayValued', ...
                'Less than two components is not supported.')
            end

            %domain will always be the same: (in polar coords)
            dom = [-pi pi 0 1]; 
            
            
            % Go pick up vectorize flag: 
            vectorize = 0; 
            for jj = 1:numel(varargin) 
                if ( strcmpi( varargin{jj}, 'vectorize' ) )
                    vectorize = 1;
                    varargin(jj) = []; 
                    break
                end
            end
            
            
            % determine if coordinate system has been provided
            % if conflicting coord systems are used, default to Cartesian.
             if isa(varargin{1}, 'diskfun') && isa(varargin{2}, 'diskfun')
                 if ((nargin < 3) && (~strcmpi(varargin{1}.coords, varargin{2}.coords)) )
                     warning('CHEBFUN:DISKFUNV:diskfunv:coords',...
                    'The two components have different coordinate settings. Now setting the diskfunv to evaluate with Cartesian coordinates.')
                    iscart = 1; 
                 else
                %if two diskfuns are given, set to match or go and get flag if present.
                 iscart = diskfun.coordsetting(varargin{:});  
                 end
             
             
             %diskfun and a function given; take diskfun setting and look
             %for flag
             elseif isa(varargin{1}, 'diskfun') && isa(varargin{2}, 'function_handle')
                 iscart = diskfun.coordsetting(varargin);
             elseif isa(varargin{1}, 'function_handle') && isa(varargin{2}, 'diskfun')
                 iscart = diskfun.coordsetting(varargin{2}, varargin{3:end});
             else %just check for flag
                 ispolar = find(strcmp(varargin,'polar'));
                    if ( any(ispolar) )
                        iscart = 0;
                    else  
                        iscart= 1; %default to Cartesian coords.
                    end
             end
            if ~iscart  %set coord sys for constructing diskfuns
                coordset = 'polar';
            else
                coordset = 'cart';
            end
            
            %get rid of flag so we can count components 
            for jj = 1:numel(varargin) 
                if ( strcmpi( varargin{jj}, 'cart' )|| strcmpi( varargin{jj}, 'polar' ) )
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
            
            %check if there was only one component
            if ( numel(fh) < 2 ) 
                error('CHEBFUN:DISKFUNV:diskfunv:arrayValued', ...
                'Less than two components is not supported.')
            end
            % Stop now if there are too many components
            if ( numel( fh ) > 2 ) 
                error('CHEBFUN:DISKFUNV:diskfunv:arrayValued', ...
                          'More than two components is not supported.')
            end 
            
            % Convert all function handles to diskfun objects: 
            for jj = 1:2
                if ( isa( varargin{jj}, 'function_handle') )
                    if ( ~vectorize )
                        newcheb = diskfun( varargin{jj}, coordset);
                    else
                        newcheb = diskfun( varargin{jj}, 'vectorize', coordset);
                    end
                    fh{jj} = newcheb;
                elseif ( isa( varargin{jj}, 'diskfun') )
                    fh{jj} = varargin{jj};
                    fh{jj}.coords = coordset;
                elseif ( isa( varargin{jj}, 'chebfun') )
                    fh{jj} = diskfunfun( varargin{jj}, coordset);
                elseif ( isa( varargin{jj}, 'double') )
                    fh{jj} = diskfun( varargin{jj}, coordset);
                end
            end
            
            % Assign to the diskfunv object: 
            F.components = fh;
            F.isTransposed = 0;
            F.domain = dom;
            F.nComponents = numel(fh);
            F.coords = coordset;
        end
    end 
    
end

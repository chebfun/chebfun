classdef spherefunv
%SPHEREFUNV   Class constructor for SPHEREFUNV objects.
% 
% SPHEREFUNV(F, G, H) constructs a SPHEREFUNV with three components from the
% function handles F, G, and H. These can also be SPHEREFUN objects or any
% other object that the SPHEREFUN constructor accepts.  Each component is
% represented as a SPHEREFUN. The domain of F, G, and H in instrinsic
% coordinates is [-pi,pi]x[0 pi].
%
% See also SPHEREFUN. 

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CLASS PROPERTIES:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    properties ( Access = public )
        components   % Array of SPHEREFUN objects.
        isTransposed % transposed?
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% STATIC METHODS:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods ( Access = public, Static = true )
        % Unit normal vector:
        n = unormal(dom);        
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CLASS CONSTRUCTOR:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods ( Access = public, Static = false )
        
        function F = spherefunv( varargin )
            % The main SPHEREFUNV constructor.
                       
            % Return an empty SPHEREFUNV:
            if ( (nargin == 0) || isempty(varargin) )
                return
            end
                       
            % This function calls the SPHEREFUN constructor once for each 
            % non-zero component because a SPHEREFUNV is just vector of 
            % SPHEREFUN objects.
            
            % If the argument is a SPHEREFUN2V, nothing to do:
            if ( isa(varargin{1}, 'spherefunv') ) 
                F = varargin{1};
                return
            end
            
            % Go pick up vectorize flag: 
            vectorize = 0; 
            for jj = 1:numel(varargin) 
                if ( strcmpi(varargin{jj}, 'vectorize') )
                    vectorize = 1;
                    varargin(jj) = []; 
                end
            end
            
            % Unwrap input arguments;
            component = 1;
            for jj = 1:numel( varargin )
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
            
            % Convert all function handles to spherefun objects: 
            for jj = 1:numel(varargin)
                if ( isa( varargin{jj}, 'function_handle') )
                    if ( ~vectorize )
                        newcheb = spherefun( varargin{jj} );
                    else
                        newcheb = spherefun( varargin{jj}, 'vectorize' );
                    end
                    fh{jj} = newcheb;
                elseif ( isa( varargin{jj}, 'spherefun') )
                    fh{jj} = varargin{jj};
                elseif ( isa( varargin{jj}, 'chebfun') )
                    fh{jj} = spherefun( varargin{jj} );
                elseif ( isa( varargin{jj}, 'double') )
                    fh{jj} = spherefun( varargin{jj} );  
                end
            end

            % Stop now if there are too many components
            if ( numel(fh) > 3 ) 
                error('CHEBFUN:SPHEREFUNV:spherefunv:arrayValued', ...
                          'More than three components is not supported.')
            end 
            
            % Stop now if there are too few components: 
            if ( numel(fh) < 3 ) 
                error('CHEBFUN:SPHEREFUNV:spherefunv:arrayValued', ...
                'Less than three components is not supported.')
            end

            % Stop now if there are no components: 
            if ( numel(fh) == 0 ) 
                error('CHEBFUN:SPHEREFUNV:spherefunv:empty', ...
                    ['The spherefunv constructor needs to be given', ...
                     'function handles or spherefun objects.'])
            end
            
            % Check the domains of all the spherfuns are the same:
            pass = zeros(numel(fh)-1,1);
            for jj = 2:numel(fh)
               pass(jj-1) = domainCheck(fh{1}, fh{jj});   
            end
            
            if ( ~all(pass) )
                error('SPHEREFUN:SPHEREFUNV:spherefunv:domainCheck', ...
                    'All spherefun objects need to have the same domain.');
            end
            
            % Assign to the spherfunv object: 
            F.components = fh;
            F.isTransposed = 0;

        end
    end 
    
end

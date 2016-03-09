%DISKFUNV   Class constructor for DISKFUNV objects.
% 
% DISKFUNV(F,G) constructs a DISKFUNV with two components from the
% function handles F and G.  These can also be DISKFUN objects or any
% other object that the DISKFUN constructor accepts.  Each component is
% represented as a DISKFUN.
%
% DISKFUNV(F,G,[A B C D]) constructs a DISKFUNV object from F, G on the domain [A B] x [C D].
%
% See also DISKFUN. 

% TO DO: add the ability to input functions in cartesian coords
% 

classdef diskfunv
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CLASS PROPERTIES:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    properties ( Access = public )
        components   % Array of DISKFUN objects.
        isTransposed % transposed?
        domain % Domain of F
        nComponents % number of components. 
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% STATIC METHODS:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods ( Access = public, Static = true )
        % Unit normal vector
        n = unormal( dom );        
    end
    
       methods ( Access = public, Static = true )
        % Unit normal vector for r (radius) and theta (angle)
        n = unit(a);        
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
            
            %domain will always be the same: (in polar coords)
            dom = [-pi pi 0 1]; 
            
            % Go and try find the domain: 
           % domain = [-pi pi 0 pi];
           % for jj = 1:numel(varargin)
              % if ( isa( varargin{jj}, 'double') && numel( varargin{jj}) == 4 ) 
              %     domain = varargin{jj}; 
               %    varargin(jj) = []; 
             %  elseif ( isa( varargin{jj}, 'spherefun') ) 
               %    domain = varargin{jj}.domain;  
              % end
            %end
            
            % Go pick up vectorize flag: 
            vectorize = 0; 
            for jj = 1:numel(varargin) 
                if ( strcmpi( varargin{jj}, 'vectorize' ) )
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
            
            % Convert all function handles to diskfun objects: 
            for jj = 1:numel(varargin)
                if ( isa( varargin{jj}, 'function_handle') )
                    if ( ~vectorize )
                        newcheb = diskfun( varargin{jj}, dom);
                    else
                        newcheb = diskfun( varargin{jj}, dom, 'vectorize');
                    end
                    fh{jj} = newcheb;
                elseif ( isa( varargin{jj}, 'diskfun') )
                    fh{jj} = varargin{jj};
                elseif ( isa( varargin{jj}, 'chebfun') )
                    fh{jj} = diskfunfun( varargin{jj}, dom);
                elseif ( isa( varargin{jj}, 'double') )
                    fh{jj} = diskfun( varargin{jj}, dom);  
                end
            end

            % Stop now if there are too many components
            if ( numel( fh ) > 3 ) 
                error('DISKFUN:DISKFUNV:diskfunv:arrayValued', ...
                          'More than three components is not supported.')
            end 
            
            % Stop now if there are too few components: 
            if ( numel( fh ) < 2 ) 
                error('DISKFUN:DISKFUNV:diskfunv:arrayValued', ...
                'Less than two components is not supported.')
            end

            % Stop now if there are no components: 
            if ( numel( fh ) == 0 ) 
                error('DISKFUN:DISKFUNV:diskfunv:empty', ...
                'The diskfun constructor needs to be given function handles or diskfun objects.')
            end
            
            % Check the domains of all the spherfuns are the same:
           % pass = zeros(numel(fh)-1,1);
           % for jj = 2:numel(fh)
            %   pass(jj-1) = domainCheck( fh{1}, fh{jj});   
          %  end
            
           % if ( ~all(pass) )
              %  error('SPHEREFUN:SPHEREFUN2V:spherefunv:domainCheck', ...
             %       'All spherefun objects need to have the same domain.');
           % end
            
            % Assign to the diskfunv object: 
            F.components = fh;
            F.isTransposed = 0;
            F.domain = dom;
            F.nComponents = numel(fh);
        end
    end 
    
end

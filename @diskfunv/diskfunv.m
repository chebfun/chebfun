%DISKFUNV   Class constructor for DISKFUNV objects.
%
%   DISKFUNV(F, G) constructs a DISKFUNV with two components from the
%   function handles F and G. It represents vector-valued functions defined
%   on the disk of radius one. These can also be DISKFUN objects. Each
%   component of a DISKFUN is represented as a DISKFUN.
%
% See also DISKFUN.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
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
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CLASS METHODS:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods

    end
    
    methods ( Access = public, Static = true )
        
        f = coeffs2diskfunv(C1, C2);    
       
        [X, Y] = coeffs2vals(U, V); 
        [U, V] = vals2coeffs(X,Y);
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
            % A. Townsend, G. Wright, H. Wilber, "Computing with functions
            % in polar and spherical geometries, II. The disk", submitted to
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
            
            %if a function handle is supplied, check for 'polar' flag
            if ( isa(varargin{1}, 'function_handle') ||...
                                     isa(varargin{2}, 'function_handle') )  
            coordset = 'cart'; 
            for jj = 1:numel(varargin)
                if ( strcmpi( varargin{jj}, 'polar' ) )
                    coordset = 'polar';
                    varargin(jj) = [];
                    break
                end
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

                    
                elseif ( isa( varargin{jj}, 'chebfun') )
                    % CHEBFUN object: 
                    fh{jj} = diskfun( varargin{jj});
                    
                elseif ( isa( varargin{jj}, 'double') )
                    % DOUBLE object:
                    fh{jj} = diskfun( varargin{jj});
                    
                end
            end
            
            % Finally, assign the two components to make a DISKFUNV object:
            F.components = fh;
            F.isTransposed = 0;
            F.domain = dom;
            F.nComponents = numel( fh );
            
        end
    end 
end
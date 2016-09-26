%CHEBFUN2V   Class constructor for CHEBFUN2V objects.
% 
% CHEBFUN2V(F,G) constructs a CHEBFUN2V with two components from the function
% handles F and G.  F and G can also be CHEBFUN2 objects or any other object
% that the CHEBFUN2 constructor accepts.  Each component is represented as a
% CHEBFUN2.
%
% CHEBFUN2V(F,G,H) constructs a CHEBFUN2V with three components from the
% function handles F, G, and H.  F, G, and H can also be CHEBFUN2 objects 
% or any other object that the CHEBFUN2 constructor accepts. 
%
% CHEBFUN2V(F,G,[A B C D]) constructs a CHEBFUN2V object from F and G 
% on the domain [A B] x [C D].
%
% CHEBFUN2V(F,G,H,[A B C D]) constructs a CHEBFUN2V object from F, G, and 
% H on the domain [A B] x [C D].
% 
% See also CHEBFUN2. 

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information. 

classdef chebfun2v
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CLASS PROPERTIES:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    properties ( Access = public )
        components   % Array of CHEBFUN2 objects.
        nComponents  % Number of components
        isTransposed % transposed?
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CLASS CONSTRUCTOR:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods ( Access = public, Static = false )
        
        function F = chebfun2v( varargin )
            % The main CHEBFUN2V constructor!
                       
            % Return an empty CHEBFUN2V:
            if ( (nargin == 0) || isempty(varargin) )
                return
            end
                       
            % This function calls the CHEBFUN2 constructor once for each 
            % non-zero component because a CHEBFUN2V is just vector of 
            % CHEBFUN2 objects.
            
            % If argument is a CHEBFUN2V, nothing to do:
            if ( isa(varargin{1}, 'chebfun2v') ) 
                F = varargin{1};
                return
            end
            
            % Go and try find the domain: 
            domain = [-1 1 -1 1]; 
            for jj = 1:numel(varargin)
               if ( isa( varargin{jj}, 'double') && numel( varargin{jj}) == 4 ) 
                   domain = varargin{jj}; 
                   varargin(jj) = []; 
               elseif ( isa( varargin{jj}, 'chebfun2') ) 
                   domain = varargin{jj}.domain;  
               end
            end
            
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
            
            % Convert all function handles to chebfun2 objects: 
            for jj = 1:numel(varargin)
                if ( isa( varargin{jj}, 'function_handle') )
                    if ( ~vectorize )
                        newcheb = chebfun2( varargin{jj}, domain);
                    else
                        newcheb = chebfun2( varargin{jj}, domain, 'vectorize');
                    end
                    fh{jj} = newcheb;
                elseif ( isa( varargin{jj}, 'chebfun2') )
                    fh{jj} = varargin{jj};
                elseif ( isa( varargin{jj}, 'double') )
                    fh{jj} = chebfun2( varargin{jj}, domain);  
                end
            end
            
            % Stop now if there are too many components
            if ( numel( fh ) > 3 ) 
                error('CHEBFUN:CHEBFUN2V:chebfun2v:arrayValued', ...
                          'More than three components is not supported.')
            end 
            
            % Stop now if there are no components: 
            if ( numel( fh ) == 0 ) 
                error('CHEBFUN:CHEBFUN2V:chebfun2v:empty', ...
                'The Chebfun2 constructor needs to be given function handles or chebfun2 objects.')
            end
            
            % Check the domains of all the chebfun2s are the same:
            pass = zeros(numel(fh)-1,1);
            for jj = 2:numel(fh)
               pass(jj-1) = domainCheck( fh{1}, fh{jj});   
            end
            
            if ( ~all(pass) )
                error('CHEBFUN:CHEBFUN2V:chebfun2v:domainCheck', ...
                    'All chebfun2 objects need to have the same domain.');
            end
            
            % Assign to the Chebfun2v object: 
            F.components = fh;
            F.nComponents = numel( fh );
            F.isTransposed = 0;

        end
    end 
    
end

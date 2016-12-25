%CHEBFUN3V   Class constructor for CHEBFUN3V objects.
% 
%   CHEBFUN3V(F, G) constructs a CHEBFUN3V with two components from the 
%   function handles F and G. F and G can also be CHEBFUN3 objects or any 
%   other object that the CHEBFUN3 constructor accepts. Each component is 
%   represented as a CHEBFUN3.
%
%   CHEBFUN3V(F, G, H) constructs a CHEBFUN3V with three components from 
%   the function handles F, G, and H.  F, G, and H can also be CHEBFUN3 
%   objects or any other objects that the CHEBFUN3 constructor accepts.
%
%   CHEBFUN3V(F, G, [A B C D E K]) constructs a CHEBFUN3V object from F and
%   G on the domain [A B] x [C D] x [E K].
%
%   CHEBFUN3V(F, G, H, [A B C D E K]) constructs a CHEBFUN3V object from F,
%   G, and H on the domain [A B] x [C D] x [E K].
% 
% See also CHEBFUN3.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

classdef chebfun3v
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CLASS PROPERTIES:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    properties ( Access = public )
        components   % Array of CHEBFUN3 objects.
        nComponents  % Number of components
        isTransposed % transposed?
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CLASS CONSTRUCTOR:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods ( Access = public, Static = false )
        
        function F = chebfun3v(varargin)
            % The main CHEBFUN3V constructor!
                       
            % Return an empty CHEBFUN3V:
            if ( (nargin == 0) || isempty(varargin) )
                return
            end
                       
            % This function calls the CHEBFUN3 constructor once for each 
            % non-zero component because a CHEBFUN3V is just vector of 
            % CHEBFUN3 objects.
            
            % If argument is a CHEBFUN3V, nothing to do:
            if ( isa(varargin{1}, 'chebfun3v') ) 
                F = varargin{1};
                return
            end
            
            % Find the domain: 
            domain = [-1 1 -1 1 -1 1]; 
            for jj = 1:numel(varargin)
               if ( isa(varargin{jj}, 'double') && numel(varargin{jj}) == 6 ) 
                   domain = varargin{jj}; 
                   varargin(jj) = []; 
               elseif ( isa( varargin{jj}, 'chebfun3') ) 
                   domain = varargin{jj}.domain;  
               end
            end
            
            % Pick up vectorize flag: 
            vectorize = 0; 
            for jj = 1:numel(varargin) 
                if ( strcmpi( varargin{jj}, 'vectorize' ) )
                    vectorize = 1;
                    varargin(jj) = []; 
                end
            end
            
            % Unwrap input arguments;
            component = 1;
            for jj = 1:numel(varargin)
                if ( iscell(varargin{jj}) ) 
                    for kk = 1:numel(varargin{jj})
                        fh{component} = varargin{jj}{kk};
                        component = component + 1; 
                    end
                else
                    fh{component} = varargin{jj};
                    component = component + 1;
                end
            end
            varargin = fh; 
            
            % Convert all function handles to chebfun3 objects: 
            for jj = 1:numel(varargin)
                if ( isa( varargin{jj}, 'function_handle') )
                    if ( ~vectorize )
                        newcheb = chebfun3(varargin{jj}, domain);
                    else
                        newcheb = chebfun3(varargin{jj}, domain, 'vectorize');
                    end
                    fh{jj} = newcheb;
                elseif ( isa(varargin{jj}, 'chebfun3') )
                    fh{jj} = varargin{jj};
                elseif ( isa(varargin{jj}, 'double') )
                    fh{jj} = chebfun3(varargin{jj}, domain);  
                end
            end
            
            % Stop if there are too many components
            if ( numel( fh ) > 3 ) 
                error('CHEBFUN:CHEBFUN3V:arrayValued', ...
                          'More than three components is not supported.')
            end 
            
            % Stop if there are no components: 
            if ( numel(fh) == 0 ) 
                error('CHEBFUN:CHEBFUN3V:empty', ...
                ['The Chebfun3 constructor needs to be given function ' ...
                     'handles or chebfun3 objects.'])
            end
            
            % Check the domains of all the chebfun3 objects are the same:
            pass = zeros(numel(fh)-1, 1);
            for jj = 2:numel(fh)
               pass(jj-1) = domainCheck(fh{1}, fh{jj});
            end
            
            if ( ~all(pass) )
                error('CHEBFUN:CHEBFUN3V:domainCheck', ...
                    'All chebfun3 objects need to have the same domain.');
            end
            
            % Assign to the Chebfun3v object: 
            F.components = fh;
            F.nComponents = numel(fh);
            F.isTransposed = 0;

        end
    end 
    
end
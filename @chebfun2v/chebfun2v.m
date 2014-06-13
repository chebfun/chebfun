% CHEBFUN2V Class constructor for CHEBFUN2V objects
% 
% CHEBFUN2V(F,G) constructs a CHEBFUN2V with two components from the function handles F
% and G.  F and G can also be CHEBFUN2 objects or any other object that the
% CHEBFUN2 constructor accepts.  Each component is represented as a CHEBFUN2. 
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

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information. 

classdef chebfun2v
    
    properties ( Access = public )
        components   % Array of CHEBFUN2 objects.
        nComponents  % Number of components
        isTransposed % transposed?
    end
    
    methods
        
        function F = chebfun2v( varargin )
            % The main CHEBFUN2V constructor!
            
            % Return an empty CHEBFUN2V:
            if ( (nargin == 0) || isempty(varargin) )
                return
            end
            
            % This function calls the CHEBFUN2 constructor once for each 
            % non-zero component because a CHEBFUN2V is just vector of 
            % CHEBFUN2 objects.
            
            op = varargin{1};
            % If argument is a CHEBFUN2V, nothing to do:
            if ( isa(op,'chebfun2v') ) 
                F = op;
                return
            end
            
            domain = [-1 1 -1 1];  % default domain
            if ( nargin >= 3 )
                if ( isa(varargin{2}, 'double') )
                    % Grab domain.
                    domain = varargin{2};
                elseif (isa( op, 'function_handle' ) &&...
                                (isa( varargin{2}, 'function_handle' ) ) )
                    % The constructor is given the components as function 
                    % handles.
                    len = numel(varargin);
                    j = 1;
                    while ( ( j < len) &&...
                                 isa( varargin{j+1}, 'function_handle' ) )
                        j = j + 1;
                    end
                    op = {op, varargin{1:j}};
                    if (j < len)
                        domain = varargin{j+1};
                    end
                    if ( j + 1 < len )
                        varargin = varargin{j+2:len};
                    else
                        varargin = [];
                    end
                elseif ( isa( op, 'chebfun2' ) &&...
                                        isa( varargin{2}, 'chebfun2' ) )
                    % The constructor is given the components as CHEBFUN2 
                    % objects.
                    domain = op.domain;
                    len = numel(varargin);
                    j = 1;
                    while ( ( j < len) && isa( varargin{j+1}, 'chebfun2' ) )
                        if ( ~domainCheck(op, varargin{j+1}) )
                            error('CHEBFUN2V:CONSTRUCTOR',...
                                                   'Domain not the same.');
                        end
                        j = j + 1;
                    end
                    op = {op, varargin{1:j}};
                    if (j < len)
                        varargin = varargin{j+1:len};
                    else
                        varargin = [];
                    end
                else
                    error('CHEBFUN2V:CONSTRUCTOR','Unrecognised syntax.');
                end
            end
             
            % Pull out domain
            for j = 1:numel( varargin )
                if ( isa( varargin{j}, 'chebfun2' ) )
                    domain = varargin{j}.domain;
                end
            end
            
            % Loop through the cell array of components.
            F.components = [];
            for j = 1:numel( varargin )
                f = chebfun2(varargin{ j }, domain, varargin );
                F.components{j} = f;
            end
            F.nComponents = numel( varargin );
            F.isTransposed = 0;
            
            if ( numel(varargin) > 3)
                error('CHEBFUN2:CONSTRUCTOR:ARRAYVALUED',...
                          'More than three components is not supported.')
            end
            
        end
    end 
    
end

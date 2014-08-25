classdef mapping
%MAPPING   Class for mapping [-1, 1] to arbitrary domains.
%
%   The MAPPING class handles the (possibly nonlinear) mapping between
%   CLASSICFUN objects and ONEFUN objects.
%
%   A MAPPING object, M, _must_ contain:
%       * M.For: The forward map from [-1, 1] to [a, b].
%       * M.Der: The derivative of M.For.
%       * M.Inv: The inverse mapping from [a, b] to [-1, 1].
%   
% See also CLASSICFUN, ONEFUN.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEVELOPER NOTE:
%   Class diagram: [<<classicFun>>] <>-- [mapping]
%                                   <>-- [oneFun]
%   Note. Currently this class onlt maps from [-1, 1] to [a, b]. This may be
%   changed in the future.
%
%   Note. We don't bother to store the target domain as this can be found by
%   evaluating the map at [-1, 1].
%
%   Note. Currently the only client is CLASSICFUN for mapping ONEFUNS. In the
%   future we might allow singular maps for dealing with singular endpoints.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    properties ( GetAccess = 'public', SetAccess = 'public' )

        % Note, we must use CAPS otherwise we cannot have 'for' as a property.
         
        % Forward map:
        For
        
        % Derivative of the map:
        Der
        
        % Inverse of the map:
        Inv

    end
    
    methods
        
        %% %%%%%%%%%%%%%%%%%%%%%%%%%%%  CONSTRUCTOR  %%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function obj = mapping(For, Der, Inv)
            % MAPPING CLASS CONSTRUCTOR.
            
            if ( nargin == 0 )
                % Return an empty MAPPING:
                return
            end
            
            if ( (nargin == 1) && isa(For, 'mapping') )
                % A MAPPING was given:
                obj = For;
                return
            end
            
            % Fill the MAPPING object. Note all fields must be given.
            obj.For = For;
            obj.Der = Der;
            obj.Inv = Inv;
            
        end
        
        function out = feval(map, x)
        %FEVAL   Evaluate the orward part of a map.
            out = feval(map.For, x);
        end
        
        function out = inv(map)
        %INV   Invert a map.
            out.For = map.Inv;
            out.Der = @(x) 1./map.Der(x);
            out.Inv = map.For;
        end
                
        % TODO: Implement this.
%         function out = compose(map1, map2)
%             % Compose two mappings.
%         end

        % TODO: Testing shows SUBSREF is expensive here. Since we don't need it,
        % (the default suffices for now) we comment it out.
%         function out = subsref(map, index) 
%             % TODO: Allow ".for", etc, in the place of ".For"?
%             idx = index(1).subs;
%             switch index(1).type
%                 case '()'
%                     if ( isa(idx{1}, 'mapping') )
%                         out = compose(map, idx{1});
%                     else
%                         out = feval(map.For, idx{:});
%                     end
%                 case '.'
%                     if ( isprop(map, idx) )
%                         out = map.(idx);
%                     elseif ( isfield(map.OtherData, idx) )
%                         out = map.OtherData.(idx);
%                     else
%                         error('CHEBFUN:MAPPING:subsref:unknown', ...
%                             'Unknown property %s.', idx);
%                     end  
%             end             
%             if ( numel(index) > 1 )
%                 out = subsref(out, index(2:end));
%             end
%         end
        
    end
    
    methods ( Static = true )
        
        function map = linear(ends)
        %LINEAR   Creates a linear map structure for BNDFUN objects.
        %   MAP = LINEAR(ENDS), where ENDS is a two-vector, returns a structure
        %   that defines a linear map. The structure MAP consists of three
        %   function handles:
        %      MAP.FOR is a function that maps [-1,1] to [ENDS(1), ENDS(2)].
        %      MAP.DER is the derivative of the map defined in MAP.FOR
        %      MAP.INV is the inverse map.
        
            a = ends(1);
            b = ends(2);
            ForHandle = @(y) b*(y + 1)/2 + a*(1 - y)/2; 
            DerHandle = @(y) (b - a)/2 + 0*y;
            InvHandle = @(x) (x - a)/(b - a) - (b - x)/(b - a);
            map = mapping(ForHandle, DerHandle, InvHandle);
            
        end
        
        function map = unbounded(ends)
        %UNBOUNDED   Creates a map structure for UNBNDFUN objects.
        %   M = UNBOUNDED(ENDS) returns a structure that defines a nonlinear map
        %   from [-1 1] to the unbounded domain [ENDS(1) ENDS(2)]. The structure
        %   MAP consists of three function handles, one string. 
        %    M.FOR is a function that maps [-1,1] to [ENDS(1) ENDS(2)]. 
        %    M.DER is the derivative of the map defined in MAP.FOR.
        %    M.INV is the inverse map of M.FOR. 
        
            
            % The domain:
            a = ends(1); 
            b = ends(2);
            
            % Fixed map parameters:
            s = 1;
            c = 0;
            
            % Deal with the different cases:
            if ( a == -inf && b == inf )
                
                ForHandle = @(y) 5*s*y./(1 - min(y.^2, 1)) + c;
                InvHandle = @(x) 2*x./(5*s + sqrt(25*s^2 + 4*x.^2));
                DerHandle = @(y) 5*s*(1 + y.^2)./(1 - y.^2).^2;
                
            elseif ( a == -inf )
                
                ForHandle = @(y) 15*s*(y - 1)./(y + 1) + b;
                InvHandle = @(x) (15*s + x - b)./(15*s - x + b);
                DerHandle = @(y) 15*s*2./(y + 1).^2;
                
            elseif ( b == inf )
                
                ForHandle = @(y) 15*s*(y + 1)./(1 - y) + a;
                InvHandle = @(x) (-15*s + x - a)./(15*s + x - a);
                DerHandle = @(y) 15*s*2./(y - 1).^2;
                
            else
                
                 error('CHEBFUN:UNBNDFUN:createMap:input', ...
                     'Error: Check input.')
                
            end
            
            map = mapping(ForHandle, DerHandle, InvHandle);
            
        end
        
    end
    
end
    

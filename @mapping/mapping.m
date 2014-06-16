classdef mapping
%MAPPING   Class for mapping [-1, 1] to arbitrary domains.
%
%   The MAPPING class handles the (possibly nonlinear) mapping between
%   CLASSICFUN objects and ONEFUN objects.
%
%   A MAPPING object, M, will contain:
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
%   evaluating the map at [-1,1].
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
        
        % Misc. data:
        OtherData = struct();
        
    end
    
    methods
        
        %% %%%%%%%%%%%%%%%%%%%%%%%%%%%  CONSTRUCTOR  %%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function obj = mapping(mapStruct)
            
            if ( nargin == 0 )
                return
            end
            
            for inProp = fieldnames(mapStruct).'
                for prop = properties(mapping()).'
                    if ( strcmpi(inProp, prop) )
                        obj.(char(prop)) = mapStruct.(char(inProp));
                        mapStruct = rmfield(mapStruct, inProp);
                    end
                end
            end
            
            if ( ~isempty(obj.OtherData) )
                obj.OtherData = mapStruct;
            end
            
        end
        
        function out = feval(map, x)
            out = feval(map.for, x);
        end
        
        function out = inv(map)
            out.For = map.Inv;
            out.Der = @(x) 1./map.Der(x);
            out.Inv = map.For;
        end
                
%         function out = compose(map1, map2)
%             % Compose two mappings.
%         end


        % TODO: Testing shows SUBSREF is expensive here. Since we don't need it,
        % (the default suffices for now) we comment it out.
%         function out = subsref(map, index) 
%             % TODO: Allow ".for", ect, in the place of ".For"?
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
        %      MAP.INV is the inverse map.
        %      MAP.DER is the derivative of the map defined in MAP.FOR
        
            a = ends(1);
            b = ends(2);
            mapStruct = struct( ...
                'For', @(y) b*(y + 1)/2 + a*(1 - y)/2, ...
                'Inv', @(x) (x - a)/(b - a) - (b - x)/(b - a), ...
                'Der', @(y) (b - a)/2 + 0*y);
            map = mapping(mapStruct);
            
        end
        
        function map = unbounded(ends)
        %UNBOUNDED   Creates a map structure for UNBNDFUN objects.
        %   M = UNBOUNDED(ENDS) returns a structure that defines a nonlinear map
        %   from [-1 1] to the unbounded domain [ENDS(1) ENDS(2)]. The structure
        %   MAP consists of three function handles, one string. M.FOR is a
        %   function that maps [-1,1] to [ENDS(1) ENDS(2)]. M.INV is the inverse
        %   map. M.DER is the derivative of the map defined in MAP.FOR.
            
            % The domain:
            a = ends(1); 
            b = ends(2);
            
            % Initialise the map structure:
            map = struct('For', [], 'Inv', [], 'Der', []);
            
            % Fixed map parameters:
            s = 1;
            c = 0;
            
            % Deal with the different cases:
            if ( a == -inf && b == inf )
                
                map.For = @(y) 5*s*y./(1 - min(y.^2, 1)) + c;
                map.Inv = @(x) 2*x./(5*s + sqrt(25*s^2 + 4*x.^2));
                map.Der = @(y) 5*s*(1 + y.^2)./(1 - y.^2).^2;
                
            elseif ( a == -inf )
                
                map.For = @(y) 15*s*(y - 1)./(y + 1) + b;
                map.Inv = @(x) (15*s + x - b)./(15*s - x + b);
                map.Der = @(y) 15*s*2./(y + 1).^2;
                
            elseif ( b == inf )
                
                map.For = @(y) 15*s*(y + 1)./(1 - y) + a;
                map.Inv = @(x) (-15*s + x - a)./(15*s + x - a);
                map.Der = @(y) 15*s*2./(y - 1).^2;
                
            else
                
                error('CHEBFUN:unbounded:input', 'Error: Check input.')
                
            end
            
        map = mapping(map);
            
        end
        
    end
    
end
    

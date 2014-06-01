classdef mapping
%MAPPING
%
%

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.
    
    properties ( GetAccess = 'public', SetAccess = 'public' )

        % Note, we must use CAPS otherwise we cannot have 'for' as a property.
         
        % Forward map:
        For
        
        % Derivative of the map:
        Der
        
        % Inverse of the map:
        Inv
        
        % Input domain:
        InDom = [-1, 1];
        
        % Output domain:
        OutDom = [-1, 1];
        
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
            out.InDom = map.OutDom;
            out.OutDom = map.InDom;
        end
        
                
%         function out = compose(map1, map2)
%             % Compose two mappings.
%         end

        
        function out = subsref(map, index) 
            
            % TODO: Allow ".for", ect, in the place of ".For"?
            
            idx = index(1).subs;
            switch index(1).type
                case '()'
                    if ( isa(idx{1}, 'mapping') )
                        out = compose(map, idx{1});
                    else
                        out = feval(map.For, idx{:});
                    end
                case '.'
                    if ( isprop(map, idx) )
                        out = map.(idx);
                    elseif ( isfield(map.OtherData, idx) )
                        out = map.OtherData.(idx);
                    else
                        error('unknown property %s', idx);
                    end
                    
            end
            
            if ( numel(index) > 1 )
                out = subsref(out, index(2:end));
            end
        end
        
    end
    
    methods ( Static = true )
        
        function map = linear(ends)
        %LINEAR   Creates a linear map structure for BNDFUN objects.
        %   MAP = LINEAR(ENDS), where ENDS is a two-vector, returns a structure that
        %   defines a linear map. The structure MAP consists of three function handles:
        %      MAP.FOR is a function that maps [-1,1] to [ENDS(1), ENDS(2)].
        %      MAP.INV is the inverse map.
        %      MAP.DER is the derivative of the map defined in MAP.FOR
        
            a = ends(1);
            b = ends(2);
            mapStruct = struct('for', @(y) b*(y + 1)/2 + a*(1 - y)/2, ...
                'inv', @(x) (x - a)/(b - a) - (b - x)/(b - a), ...
                'der', @(y) (b - a)/2 + 0*y);
            map = mapping(mapStruct);
            
        end
        
        function map = unbounded(ends)
        %UNBOUNDED   Creates a map structure for UNBNDFUN objects.
        %   M = UNBOUNDED(ENDS) returns a structure that defines a nonlinear map
        %   from [-1 1] to the unbounded domain [ENDS(1) ENDS(2)].
        %   The structure MAP consists of three function handles, one string.
        %   M.FOR is a function that maps [-1,1] to [ENDS(1) ENDS(2)].
        %   M.INV is the inverse map.
        %   M.DER is the derivative of the map defined in MAP.FOR.
            
            % The domain:
            a = ends(1); 
            b = ends(2);
            
            % Initialise the map structure:
            map = struct('for', [], 'inv', [], 'der', [], 'forDerExps', []);
            
            % Fixed map parameters:
            s = 1;
            c = 0;
            
            % Deal with the different cases:
            if ( a == -inf && b == inf )
                
                map.for = @(y) 5*s*y./(1 - min(y.^2, 1)) + c;
                map.inv = @(x) 2*x./(5*s + sqrt(25*s^2 + 4*x.^2));
                map.der = @(y) 5*s*(1 + y.^2)./(1 - y.^2).^2;
                map.forDerExps = [-2 -2];
                
            elseif ( a == -inf )
                
                map.for = @(y) 15*s*(y - 1)./(y + 1) + b;
                map.inv = @(x) (15*s + x - b)./(15*s - x + b);
                map.der = @(y) 15*s*2./(y + 1).^2;
                map.forDerExps = [-2 0];
                
            elseif ( b == inf )
                
                map.for = @(y) 15*s*(y + 1)./(1 - y) + a;
                map.inv = @(x) (-15*s + x - a)./(15*s + x - a);
                map.der = @(y) 15*s*2./(y - 1).^2;
                map.forDerExps = [0 -2];
                
            else
                
                error('CHEBFUN:unbounded:input', 'Error: Check input')
                
            end
            
        map = mapping(map);
            
        end
        
    end
    
end
    

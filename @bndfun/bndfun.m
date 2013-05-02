classdef bndfun < fun
%BNDFUN  BNDFUN class for representing global functions on bounded [a, b].

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BNDFUN Class Description:
%
% The FUN class is an abstract class for representations of functions on the
% finite interval [a, b]. It achieves this by taking a onefun on [-1, 1] and
% applying a linear mapping.
%
% Class diagram: [<<fun>>] <-- [<<bndfun>>] <>-- [<<onefun>>]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    methods
        
        function obj = bndfun(op, domain, hscale, vscale, pref)
            
            % Construct an empty bndfun
            if ( nargin == 0 || isempty(op) )
                return
            end
            
            % Obtain preferences if none given:
            if ( nargin < 5 )
                pref = bndfun.pref;
            end
            
            % Use default domain if none given:
            if ( nargin < 2 || isempty(domain) )
                domain = pref.bndfun.domain;
            end
            
            % Check domain
            if ( any(isinf(domain)) )
                error('CHEBFUN:BNDFUN:UNBND',...
                    'Should not encounter unbounded domain in bndfun class.');
            end
            
            % Define scales if none given.
            if ( nargin < 3 )
                hscale = norm(domain,inf); 
            end
            if ( nargin < 4 )
                vscale = 0; 
            end

            linmap = bndfun.createMap(domain);
            % Include linear mapping from [-1,1] to [a,b] in the op:
            if ( isa(op, 'function_handle') && ~all(domain == [-1 1]) && ~isnumeric(op) )
                op = @(x) op(linmap.for(x));
            end
            
            % Call the OneFun constructor:
            pref = onefun.pref(pref, pref.bndfun);
            obj.onefun = onefun.constructor(op, vscale, hscale, pref);
            
            % Add the domain and mapping:
            obj.domain = domain;
            obj.mapping = linmap;
            
        end

    end
    
    % Static methods implimented by BNDFUN class.
    methods ( Static = true ) 
        
        % Retrieve and modify preferences for this class.
        prefs = pref(varargin)
        
        % Linear map from [-1, 1] to the domain of the BNDFUN.
        m = createMap(domain);  
    end
    
    
end
   
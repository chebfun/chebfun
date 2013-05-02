classdef unbndfun < fun
%UNBNDFUN UNBNDFUN class for representing global functions on unbounded domains.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% UNBNDFUN Class Description:
%
% The UNBNDFUN class is an abstract class for representations of global
% functions on the infinite interval [-inf, b], [a, inf], or [-inf, inf]. It
% achieves this by taking a onefun on [-1, 1] and applying a nonlinear mapping.
%
% Class diagram: [<<fun>>] <-- [unbndfun] <>-- [<<onefun>>]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    methods
        
        function obj = unbndfun(op, domain, hscale, vscale, pref)
            
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
            if ( ~any(isinf(domain)) )
                error('CHEBFUN:UNBNDFUN:BoundedDomain',...
                    'Should not encounter bounded domain in bndfun class.');
            end
            
            % Define scales if none given.
            if ( nargin < 3 )
                % The hscale of an unbounded domain is alway 2.
                % [TODO]: Why!?
                hscale = 2;
            end
            if ( nargin < 4 )
                vscale = 0; 
            end

            m = obj.createMap(domain);
            % Include nonlinear mapping from [-1,1] to [a,b] in the op:
            if ( isa(op, 'function_handle') && ~all(domain == [-1 1]) ...
                                            && ~isnumeric(op) )
                op = @(x) op(m.for(x));
            elseif ( isnumeric(op) )
                % [TODO]: Does this make sense for an UNBNDFUN?
            end
            
            % Call the ONEFUN constructor:
            pref = onefun.pref(pref, pref.bndfun);
            obj.onefun = onefun.constructor(op, vscale, hscale, pref);
            
            % Add the domain and mapping:
            obj.domain = domain;
            obj.mapping = m;
            
        end

    end
    
    % Static methods implimented by UNBNDFUN class.
    methods ( Static = true ) 
        
        % Retrieve and modify preferences for this class.
        prefs = pref(varargin)
        
        % Make method
        f = make(varargin)
        
        % Noninear map from [-1, 1] to the domain of the UNBNDFUN.
        m = createMap(domain);  
        
    end
    
    
end
   
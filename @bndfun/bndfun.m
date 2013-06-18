classdef bndfun < fun
%BNDFUN   Represent globally smooth functions on a bounded interval [a, b].

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BNDFUN Class Description:
%
% The BNDFUN class is a class for representations of functions on the finite
% interval [a, b]. It achieves this by taking a ONEFUN on [-1, 1] and applying
% a linear mapping.
%
% Note that all binary BNDFUN operators (methods which can take two BNDFUN
% arguments) assume that the domains of the BNDFUN objects agree. The methods
% will not issue warnings if this condition is violated, but the results will
% not be meaningful.
%
% Class diagram: [<<fun>>] <-- [<<bndfun>>] <>-- [<<onefun>>]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %% CLASS CONSTRUCTOR:
    methods
        
        function obj = bndfun(op, domain, vscale, hscale, pref)
            
            % Construct an empty bndfun:
            if ( nargin == 0 || isempty(op) )
                return
            end
            
            % Obtain preferences if none given:
            if ( nargin < 5 || isempty(pref))
                pref = bndfun.pref;
            else
                pref = bndfun.pref(pref);
            end
            
            % Use default domain if none given:
            if ( nargin < 2 || isempty(domain) )
                domain = pref.bndfun.domain;
            end
            
            % Check domain:
            if ( any(isinf(domain)) )
                error('CHEBFUN:BNDFUN:UNBND',...
                    'Should not encounter unbounded domain in bndfun class.');
            end
            
            % Define scales if none given:
            if ( nargin < 3 || isempty(vscale) )
                vscale = norm(domain,inf); 
            end
            if ( nargin < 4 || isempty(hscale) )
                hscale = 0; 
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
    
    %% STATIC METHODS IMPLEMENTED BY BNDFUN CLASS.
    methods ( Static = true ) 
        
        % Retrieve and modify preferences for this class.
        prefs = pref(varargin)
        
        % Linear map from [-1, 1] to the domain of the BNDFUN.
        m = createMap(domain);
        
        % Make a BNDFUN (constructor shortcut):
        f = make(varargin);
        
    end
    
    %% METHODS IMPLEMENTED BY THIS CLASS.
    methods
        
        % Compose a BNDFUN with an operator or another BNDFUN
        f = compose(f, op, g, pref)
        
        % Indefinite integral of a BNDFUN.
        f = cumsum(f, m, pref)
        
        % Derivative of a BNDFUN.
        f = diff(f, k, dim)
       
        % Change of domains of BNDFUN via linear change of variables
        f = changeMap(f,newdom)
        
        % Flip/reverse a BNDFUN object.
        f = flipud(f)
        
        % QR factorisation of an array-valued BNDFUN.
        [f, R, E] = qr(f, flag)

        % Restrict a BNDFUN to a subinterval.
        f = restrict(f, s)
        
        % Definite integral of a BNDFUN on the interval [a,b].
        out = sum(f, dim)
    end
end
   

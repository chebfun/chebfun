classdef fun % (Abstract)
%FUN  Abstract FUN class for representing global functions on [a, b].

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FUN Class Description:
%
% The FUN class is an abstract class for representations of functions on the
% interval [a, b]. It acheives this my taking a onefun on [-1, 1] and applying
% a mapping.
%
% The current instances of FUNs are BNDFUNS and UNBNDFUNS. The former are used
% to represent functions on bounded domains, whereas the latter are able to
% represent some functions on unboujnded domains.
%
% Class diagram: [chebfun] <>-- [<<FUN>>] <>-- [<<onefun>>]
%                                 ^   ^
%                                /     \
%                          [bndfun]   [unbndfun]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    properties (Access = public)
        domain
        mapping
        onefun
    end
    
    methods (Static)
        
        function obj = constructor(op, domain, hscale, vscale, pref)
            
            % Construct an empty fun:
            if ( nargin == 0 )
                op = [];
            end
            
            % Obtain preferences if none given:
            if ( nargin < 5 )
                pref = fun.pref;
            else
                pref = fun.pref(pref);
            end
           
            % Get domain if none given:
            if ( nargin < 2 )
                domain = pref.fun.domain;
            end

            % Get scales if none given:
            if ( nargin < 3 || isstruct(hscale) )
                if ( nargin > 2 && isstruct(hscale) )
                    pref = hscale; 
                end
                hscale = norm(domain, inf); 
                if ( isinf(hscale) )
                    hscale = 1; 
                end
            end
            if ( nargin < 4 )
                vscale = 0;
            end

            % Call constructor depending on domain:
            if ( ~any(isinf(domain)) )
                pref = bndfun.pref(pref, pref.fun);
                obj = bndfun(op, domain, hscale, vscale, pref);
            else
                pref = unbndfun.pref(pref, pref.fun);
                obj = unbndfun(op, domain, hscale, vscale, pref);
            end
            
        end

    end
    
    % Static methods implimented by FUN class.
    methods (Static = true)
        
        % Retrieve and modify preferences for this class.
        prefs = pref(varargin);
        
        % Linear map from [-1, 1] to bounded domain.
        m = linear(domain);  
        
        % Edge detector.
        [edge, vscale] = detectedge(op, domain, hscale, vscale, pref, d)

    end

end

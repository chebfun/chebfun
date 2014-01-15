classdef fun % (Abstract) 
%FUN   Approximate functions on arbitrary domains.
%   Abstract (interface) class for approximating functions on the arbitrary 
%   intervals.
%
% See also DELTAFUN, CLASSICFUN.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FUN Class Description:
%  [TODO]
%
% Class diagram: [<<CHEBFUN>>] <>-- [<<FUN>>] <----[<<classicfun>>]
%                                             <----[    deltafun  ]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Copyright 2013 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

    methods (Static)
        function obj = constructor(op, domain, vscale, hscale, pref)            
            
            if ( nargin == 0  )
                op = [];
            end
            
            % Obtain preferences if none given:
            if ( nargin < 5 )
                pref = chebpref();
            else
                pref = chebpref(pref);
            end
            
            % Check if delta functions are required:
            if ( pref.enableDeltaFunctions )                
                % Generalized function mode; call DELTAFUN constructor:                
                % Then op is a classicfun, vscale and hscale are magnitude 
                % and location of delta functions. domain is a spurious argument.
                deltaMag = vscale;
                deltaLoc = hscale;
                %[TODO]: pass preferences as well.
                obj = deltafun(op, deltaMag, deltaLoc);
            else             
                % Get domain if none given:
                if ( nargin < 2 || isempty(domain) )
                    domain = pref.domain;
                end

                % Get vscale if none given:
                if ( nargin < 3 || isstruct(vscale) )
                    vscale = 0;
                end

                % Get hscale if none given:
                if ( nargin < 4 || isempty(vscale) )
                    hscale = norm(domain, inf);
                end

                % [TODO]: Explain this. Only becomes relevant with UNBNDFUN
                if ( isinf(hscale) )
                    hscale = 1;
                end

                % Call the relevent constructor:
                if ( isa(op, 'fun') )
                    % OP is already a ONEFUN!
                    obj = op;               
                else
                    % STANDARD mode; call SMOOTHFUN constructor:
                    obj = classicfun.constructor(op, domain, vscale, hscale, pref);

                end
            end
        
        end
        
    end
    
    %% ABSTRACT (NON-STATIC) METHODS REQUIRED BY FUN CLASS.
    methods ( Abstract = true )

    end

    %% ABSTRACT STATIC METHODS REQUIRED BY ONEFUN CLASS.
    methods ( Abstract = true, Static = true )

    end
    
    %% Methods implemented by ONEFUN class.
    methods 
        
    end
    
    %% Static methods implemented by ONEFUN class.
    methods ( Static = true ) 
        
    end
    
end

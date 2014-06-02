classdef fun % (Abstract) 
%FUN   Approximate functions on arbitrary domains.
%   Abstract (interface) class for approximating functions on the arbitrary 
%   intervals.
%
% See also DELTAFUN, CLASSICFUN.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FUN Class Description:
%  [TODO]
%
% Class diagram: [<<CHEBFUN>>] <>-- [<<FUN>>] <----[<<classicfun>>]
%                                             <----[    deltafun  ]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Copyright 2014 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

    methods (Static)
        function obj = constructor(op, data, pref)
            
            % We can't return an empty FUN, so pass an empty OP down.
            if ( nargin == 0  )
                op = [];
            end

            % Parse inputs.
            if ( (nargin < 2) || isempty(data) )
                    data = struct();
            end

            if ( (nargin < 3) || isempty(pref) )
                pref = chebfunpref();
            else
                pref = chebfunpref(pref);
            end

            [domain, hscale] = parseDataInputs(data, pref);

            % [TODO]: Explain this. Only becomes relevant with UNBNDFUN
            if ( isinf(hscale) )
                data.hscale = 1;
            end
            
            % Check if delta functions are required:
            if ( pref.enableDeltaFunctions && isfield(data, 'deltaLoc') )
                % Generalized function mode; call DELTAFUN constructor:                
                obj = deltafun(op, data, pref);
            elseif ( isa(op, 'fun') )             
                % OP is already a ONEFUN!
                obj = op;               
            else
                % STANDARD mode; call SMOOTHFUN constructor:
                obj = classicfun.constructor(op, data, pref);
            end
        end
        
    end
    
    %% ABSTRACT (NON-STATIC) METHODS REQUIRED BY FUN CLASS.
    methods ( Abstract = true )

    end

    %% ABSTRACT STATIC METHODS REQUIRED BY THIS CLASS.
    methods ( Abstract = true, Static = true )

        % Map from [-1, 1] to the domain of the FUN.
        m = createMap(domain);  
        
        % Make a FUN. (Constructor shortcut)
        f = make(varargin);
    end
    
    methods ( Abstract = true, Static = false )
        
    end
    
end

function [domain, hscale] = parseDataInputs(data, pref)

domain = getDataInput(data, 'domain',  pref.domain);
hscale = getDataInput(data, 'hscale',  norm(domain, inf));

end

function val = getDataInput(data, field, defaultVal)

if ( isfield(data, field) && ~isempty(data.(field)) )
    val = data.(field);
else
    val = defaultVal;
end

end

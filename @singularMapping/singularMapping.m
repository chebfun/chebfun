classdef singularMapping < mapping
%TODO: Document    
    
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%  PROPERTIES  %%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    properties ( GetAccess = 'public', SetAccess = 'public' )
        
        % Parameters:
        params
        
    end
    
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%  METHODS  %%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    methods
        
        function obj = singularMapping(For, Der, Inv, params)
            
            obj@mapping(For, Der, Inv);
            obj.params = params;
            
        end
        
        function out = isequal(f,g)
            
            if ( ~isa(f, 'singularMapping') || ~isa(g, 'singularMapping') )
                out = false;
            elseif ( all(f.params == g.params) )
                out = true;
            else
                out = false;
            end
            
        end
        
        function out = islinear(f)
            out = false;
        end  
        
    end
    
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%  STATIC METHODS  %%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    methods ( Static = true )
        
        function m = singmap(ends, pow)
            % TODO: Document (taken from chebfun v4)
            
            default_pow = .25;
            if ( nargin == 1 )
                pow = [default_pow, default_pow];
            end
            
            pow(~pow) = 1;
            
            pos = 0;
            for k = 1:2
                if ( pow(k)~=1 && pow(k)~=0 )
                    pos = pos+(-1)^k;
                elseif ( isnan(pow(k)) )
                    pow(k) = default_pow;
                end
            end
            
            L = mapping.linear(ends);
            
            % The map is linear
            if ( ~pos && all(pow == 1) )
                m = L;
                return
            end
            
            % Can only do .25 powers at boths ends (or maybe .5?)
            if ( ~pos && (~all(pow==.25) && ~all(pow==.5)) )
                pow = [.25, .25];
                warning('FUN:sing:bothends',['Singmaps at boths ends may only have ', ...
                    'parameters 0.25 or 0.5']);
                
            end
            
            powi = 1./pow;
            
            switch pos
                
                case -1 % Left point singularity
                    powi = powi(1);
                    For = @(y) L.For( 2*( .5*(y+1) ).^powi - 1 );
                    %         m.der = @(y) L.Der(1) * 2 * powi * ( .5*(y+1) ).^(powi-1);
                    Inv = @(x) 2*( .5*(L.Inv(x)+1) ).^pow(1) - 1;
                    Der = @(y) L.Der(1) * 1 * powi * ( .5*(y+1) ).^(powi-1);
                case 1 % Right point singularity
                    powi = powi(2);
                    For = @(y) L.For( 1 - 2*( .5*(1-y) ).^powi);
                    %         m.der = @(y) -L.Der(1) * 2 * powi * ( .5*(1-y) ).^(powi-1);
                    Inv = @(x) 1 - 2*( .5*(1-L.Inv(x)) ).^pow(2);
                    Der = @(y) L.Der(1) * 1 * powi * ( .5*(1-y) ).^(powi-1);
                case 0 % Both points sigularities
                    if ( all(pow(1) == .5) )
                        For = @(y) L.For(sin(pi/2*y));
                        Inv = @(x) 2/pi*asin(L.Inv(x));
                        Der = @(y) L.Der(1)*(pi/2)*cos(pi/2*y);
                    else
                        For = @(y) L.For(sin(pi/2*sin(pi/2*y)));
                        Inv = @(x) 2/pi*asin(2/pi*asin(L.Inv(x)));
                        Der = @(y) L.Der(1)*(1/4)*cos((1/2)*pi*sin((1/2)*pi*y)).*pi^2.*cos((1/2)*pi*y);
                    end
            end
            
            m = singularMapping(For, Der, Inv, pow);
            
        end
        
        
    end
    
end
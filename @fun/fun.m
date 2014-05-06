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
        function obj = constructor(op, domain, vscale, hscale, pref)            
            
            if ( nargin == 0  )
                op = [];
            end
            
            % Obtain preferences if none given:
            if ( nargin < 5 )
                pref = chebfunpref();
            else
                pref = chebfunpref(pref);
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

    %% ABSTRACT STATIC METHODS REQUIRED BY THIS CLASS.
    methods (Abstract = true, Static = true)

        % Map from [-1, 1] to the domain of the FUN.
        m = createMap(domain);  
        
        % Make a FUN. (Constructor shortcut)
        f = make(varargin);
    end
    
    %% METHODS IMPLEMENTED BY THIS CLASS.
    methods
        
        % Absolute value of a FUN. (f should have no zeros in its domain)
        f = abs(f, pref)

        % FUN logical AND.
        h = and(f, g)

        % True if any element of a FUN is a nonzero number, ignoring NaN.
        a = any(f, dim)

        % Plot (semilogy) the Chebyshev coefficients of a FUN object, if it is
        % based on Chebyshev technology.
        h = chebpolyplot(f, varargin)

        % Complex conjugate of a FUN.
        f = conj(f)
        
        % FUN objects are not transposable.
        f = ctranspose(f)
        
        % Extract information for DISPLAY.
        info = dispData(f)
        
        % Extract boundary roots and represent them by an appropriate ONEFUN.
        f = extractBoundaryRoots(f)

        % Extract columns of an array-valued FUN object.
        f = extractColumns(f, columnIndex);

        % Round a FUN towards zero.
        g = fix(f);
        
        % Round a FUN towards minus infinity.
        g = floor(f);

        % Flip columns of an array-valued FUN object.
        f = fliplr(f)
        
        % Get properties of a FUN.
        out = get(f, prop);
        
        % Imaginary part of a FUN.
        f = imag(f)

        % Test if a FUN object manages delta functions.
        out = isdelta(f)

        % True for an empty FUN.
        out = isempty(f)

        % Test if FUN objects are equal.
        out = isequal(f, g)

        % Test if a FUN is bounded.
        out = isfinite(f)

        % Test if a FUN is unbounded.
        out = isinf(f)

        % Test if a FUN has any NaN values.
        out = isnan(f)

        % True for real FUN.
        out = isreal(f)
        
        % Test if a FUN object is built upon SINGFUN.
        out = issing(f)
        
        % Test if a FUN object is built upon SMOOTHFUN.
        out = issmooth(f)
        
        % Test if a FUN object is defined on an unbounded domain.
        out = isunbnd(f)

        % True for zero FUN objects.
        out = iszero(f)
        
        % Return Legendre coefficients of a FUN object.
        c_leg = legpoly(f, n)
        
        % Length of a FUN.
        len = length(f)

        % FUN logical.
        f = logical(f)

        % Convert an array-valued FUN into a cell array of FUN objects.
        g = mat2cell(f, M, N)

        % Global maximum of a FUN on [a,b].
        [maxVal, maxPos] = max(f)

        % Global minimum of a FUN on [a,b].
        [minVal, minPos] = min(f)

        % Global minimum and maximum on [a,b].
        [vals, pos] = minandmax(f)

        % Subtraction of two FUN objects.
        f = minus(f, g)

        % Multiplication of FUN objects.
        f = mtimes(f, c)
        
        % Estimate the Inf-norm of a FUN object.
        out = normest(f);
        
        % FUN logical NOT.
        f = not(f)

        % FUN logical OR.
        h = or(f, g)

        % Basic linear plot for FUN objects.
        varargout = plot(f, varargin)
    end
    
end

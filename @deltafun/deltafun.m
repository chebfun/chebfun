classdef deltafun
    %DELTAFUN   Class for distributions based on Dirac-delta functions on arbitrary
    %intervals.
    %
    %   Class for approximating generalized functions on the interval [a, b]
    %   using a chebfun part with no distributions and a singular part containing
    %   delta functions.
    %
    %   DELTAFUN class description
    %   [TODO]:
    %
    %   [TODO]: Calling Sequence
    %
    % See also PREF
    
    % Copyright 2013 by The University of Oxford and The Chebfun Developers.
    % See http://www.chebfun.org/ for Chebfun information.
    
    %% Properties of SINGFUN objects
    properties ( Access = public )
        % Smooth part of the representation.
        funPart     % (smoothfun)
        
        % Delta functions' structure
        delta = struct( 'magnitude', [], 'location', [], ...
                        'isReal', [], 'isImag', [], 'isConj', [] );
        % Domain
        domain
        
        % isTransposed flag
        isTransposed
    end
    
    %% CLASS CONSTRUCTOR:
    methods ( Static = true )
        function obj = deltafun(magnitude, location, funPart, domain, pref)
            %%
            % Check for preferences in the very beginning.
            if ( (nargin < 5) || isempty(pref) )
                % Determine preferences if not given.
                pref = deltafun.pref;
            else
                % Merge if some preferences are given.
                pref = deltafun.pref(pref);
            end
            
            
            %% Cases based on the number of arguments
            % Case 0: No input arguments, return an empty object.
            if ( nargin == 0 )
                obj.funPart = [];
                obj.delta = [];
                obj.domain = [];
                obj.isTransposed = [];
                return
            end
            
            %%
            % Case 1: One input argument.
            if ( nargin == 1 )
                error('CHEBFUN:DELTAFUN:nargin', 'At least two arguments should be given.');
            end
            %%
            % Case 2: Two input arguments.
            if ( nargin == 2 )
                % Decide domain:
                % [TODO]: At the moment we just define the domain as
                a = min(location);
                b = max(location);
                if( a == b ) % we should do something better here.
                    a = a - 1;
                    b = b + 1;
                end
                domain = [a - (b-a)/2, b + (b-a)/2];
                % Function part is zero:
                funPart = chebfun(0, domain);
                
            end
            
            %%
            % Case 3: Three or more input arguments.
            if ( nargin == 3 )
                % Domain not given, use the default domain used by the funPart:
                domain = funPart.domain;
            end
            
            %%
            if ( nargin >= 4 )
                % Do nothing here, all checks on arguments are done below.
            end
            %%
            % Various checks on argument compatibilities
            
            if ( size(magnitude, 2) ~= length(location) )
                error('CHEBFUN:DELTAFUN:dim', 'Magnitude should have the same number of columns as locations' );
            end
            
            if ( min(size(location)) > 1 )
                error('CHEBFUN:DELTAFUN:dim', 'Magnitude and location should each be a vector');
            end
            
            % There should be no duplicates in location:
            if ( numel(location) ~= numel(unique(location)) )
                error('CHEBFUN:DELTAFUN:duplication', 'No duplicates are allowed in location.');
            end
            
            % Make sure location is a row vector:
            location = location(:).';
            
            % Locations of delta functions should be within the domain:
            % NOTE: In fact they should be strictly in the interior of the
            % domain to make sense.
            if( any(size(domain) ~= [1 2]) )
                error('CHEBFUN:DELTAFUN:domain', 'domain must be a 1 x 2 row vector.' );
            end
            
            if( max(location) > domain(2) || min(location) < domain(1)  )
                error('CHEBFUN:DELTAFUN:domain', 'Location of a delta fun is outside the domain');
            end
            
            % Domains of deltaFun and the chebfun part should overlap:
            if( domain ~= funPart.domain )                
                error('CHEBFUN:DELTAFUN:domain', 'Domain of deltaFun should be the same as its funPart.');                
            end
                
            % Now that we have checked all the arguments, copy them in the
            % current object:
            obj.delta.magnitude = magnitude;
            obj.delta.location  = location;
            obj.delta.isImag    = false * location;
            obj.delta.isReal    = false * location;
            obj.delta.isConj    = false * location;
            
            obj.domain       = domain;
            obj.funPart      = funPart;
            obj.isTransposed = 0;
        end
    end
    
    %%
    
    
    %% METHODS IMPLEMENTED BY THIS CLASS.
    methods
        
        % Complex conjugate of a DELTAFUN.
        f = conj(f)
        
        % DELTAFUN obects are not transposable.
        f = ctranspose(f)
        
        % Indefinite integral of a DELTAFUN.
        f = cumsum(f, m, pref)
        
        % Derivative of a DELTAFUN.
        f = diff(f, k)
        
        % Dirac delta function. 
        d = dirac(f, k)
        
        % Evaluate a DELTAFUN.
        y = feval(f, x)
        
        % Flip columns of an array-valued DELTAFUN object.
        f = fliplr(f)
        
        % Flip/reverse a DELTAFUN object.
        f = flipud(f)
        
        % Imaginary part of a DELTAFUN.
        f = imag(f)
        
        % Innerproduct, equivalent to action of a distribution 
        % on a chebfun
        out = ip(f,g)
        
        % True for an empty DELTAFUN.
        out = isempty(f)
        
        % Test if DELTAFUN objects are equal.
        out = isequal(f, g)
        
        % Test if a SINGFUN is bounded.
        out = isfinite(f)
        
        % Test if a SINGFUN is unbounded.
        out = isinf(f)
        
        % Test if a SINGFUN has any NaN values.
        out = isnan(f)
        
        % True for real SINGFUN.
        out = isreal(f)
        
        % True if the DELTAFUN object has no delta functions       
        out = issmooth(f)
        
        % Ture if the DELTAFUN object is zero
        out = iszero(f)
        
        % Length of a SINGFUN.
        len = length(f)
        
        % Convert a array-valued SINGFUN into an ARRAY of SINGFUN objects.
        g = mat2cell(f, M, N)
        
        % Global maximum of a SINGFUN on [-1,1].
        [maxVal, maxPos] = max(f)
        
        % Global minimum of a SINGFUN on [-1,1].
        [minVal, minPos] = min(f)
        
        % Global minimum and maximum of a SINGFUN on [-1,1].
        [vals, pos] = minandmax(f)
        
        % Subtraction of two SINGFUN objects.
        f = minus(f, g)
        
        % Left matrix divide for SINGFUN objects.
        X = mldivide(A, B)
        
        % Right matrix divide for a SINGFUN.
        X = mrdivide(B, A)
        
        % Multiplication of SINGFUN objects.
        f = mtimes(f, c)
        
        % Basic linear plot for SINGFUN objects.
        varargout = plot(f, varargin)
        
        % Obtain data used for plotting a SINGFUN object.
        data = plotData(f)
        
        % Addition of two SINGFUN objects.
        f = plus(f, g)
        
        % Return the points used by the smooth part of a SINGFUN.
        out = points(f)
        
        % Dividing two SINGFUNs
        f = rdivide(f, g)
        
        % Real part of a SINGFUN.
        f = real(f)
        
        % Restrict a SINGFUN to a subinterval.
        f = restrict(f, s)
        
        % Roots of a SINGFUN in the interval [-1,1].
        out = roots(f, varargin)
        
        % Size of a SINGFUN.
        [siz1, siz2] = size(f, varargin)
        
        % Simplify a DELTAFUN
        f = simplify(f)
        
        % Definite integral of a SINGFUN on the interval [-1,1].
        out = sum(f, dim)
        
        % SINGFUN multiplication.
        f = times(f, g)
        
        % DELTAFUN objects are not transposable.
        f = transpose(f)
        
        % Unary minus of a DELTAFUN.
        f = uminus(f)
        
        % Unary plus of a DELTAFUN.
        f = uplus(f)                
    end
    
    %% STATIC METHODS IMPLEMENTED BY THIS CLASS.
    methods ( Static = true )
        % smooth fun constructor
        s = constructFunPart( op, pref)
        
        % remove zero columns
        [A, v] = cleanColumns(A, v);
        
        % remove zero trailing rows
        A = cleanRows(A);
        
        % Merge columns of a matrix based on duplicate values in v.
        [A, v] = mergeColumns(A, v)
        % Retrieve and modify preferences for this class.
        prefs = pref(varargin)
        
        % Costruct a zero DELTAFUN
        s = zeroDeltaFun(domain)
    end
    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Functions implemented in this file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

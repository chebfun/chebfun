classdef deltafun    
%DELTAFUN   Class for distributions based on Dirac-delta functions.
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
        
        % Delta functions.
        deltaMag       % (vector of magnitudes)
        
        % Delta Locaions.
        deltaLoc
        
        % Domain
        domain
        
        % Order of the derivative
        diffOrder    % (1x1 double)
        
        % If the imaginary part only is needed 
        isImag       % (1x1 logical)
        
        % If the imaginary part only is needed
        isReal       % (1x1 logical)
        
        % If the conjugate is needed
        isConj       % (1x1 logical)
        
        
     end
    
    %% CLASS CONSTRUCTOR:
    methods ( Static = true )
        function obj = deltafun(magnitude, location, funPart, pref) 
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
                obj.deltaMag = [];
                obj.deltaLoc = [];
                obj.domain = [];
                obj.diffOrder = [];
                obj.isImag = [];
                obj.isReal = [];
                obj.isConj = [];
                return
            end
            
            %%
            % Case 1: One input argument.
            if ( nargin == 1 )
                error('ahhh, one argument for deltafun, what?');                
            end
            %%
            % Case 2: Two input arguments.
            if ( nargin == 2 )
                
            end
                
            %%
            % Case 3: Three or more input arguments.
            if ( nargin == 3 )
                % Domain not given, use the default domain.
                domain = [-1, 1];
            end
            
            if ( length(magnitude) ~= length(location) )
                error('CHEBFUN:DELTAFUN:dim', 'Magnitude and location should be vectors of the same size.');
            end
                
            if ( min(size(magnitude)) > 1 || min(size(location)) > 1 )
                error('CHEBFUN:DELTAFUN:dim', 'Magnitude and location should each be a vector');
            end
            
            % There should be no duplicates.
            if ( numel(location) ~= numel(unique(location)) )
                error('CHEBFUN:DELTAFUN:duplication', 'No duplicates are allowed in location.');
            end
            
            % Make sure magnitude and location are row vectors.
            magnitude = magnitude(:).';
            location = location(:).';
            
            obj.deltaMag = magnitude;
            obj.deltaLoc = location;
            if( max(abs(location)) > 1 )
                error('CHEBFUN:DELTAFUN:domain', 'Domain not provided.');
            else
                obj.domain   = [-1, 1];
            end
            obj.diffOrder = 0*magnitude;
            obj.isImag = false * location;
            obj.isReal = false * location;
            obj.isConj = false * location;                        
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
        
        % Evaluate a DELTAFUN.
        y = feval(f, x)

        % Flip columns of an array-valued DELTAFUN object.
        f = fliplr(f)
        
        % Flip/reverse a DELTAFUN object.
        f = flipud(f)
        
        % Imaginary part of a DELTAFUN.
        f = imag(f)
     
        % True for an empty DELTAFUN.
        out = isempty(f)

        % Test if SINGFUN objects are equal.
        out = isequal(f, g)

        % Test if a SINGFUN is bounded.
        out = isfinite(f)

        % Test if a SINGFUN is unbounded.
        out = isinf(f)

        % Test if a SINGFUN has any NaN values.
        out = isnan(f)

        % True for real SINGFUN.
        out = isreal(f)
        
        % True for zero SINGFUN objects
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

        % Definite integral of a SINGFUN on the interval [-1,1].
        out = sum(f, dim)

        % SINGFUN multiplication.
        f = times(f, g)
        
        % SINGFUN objects are not transposable.
        f = transpose(f)

        % Unary minus of a SINGFUN.
        f = uminus(f)

        % Unary plus of a SINGFUN.
        f = uplus(f)
                
    end

    %% STATIC METHODS IMPLEMENTED BY THIS CLASS.
    methods ( Static = true )                                
        % smooth fun constructor
        s = constructFunPart( op, pref)
                
        % Retrieve and modify preferences for this class.
        prefs = pref(varargin)
        
        % Costruct a zero SINGFUN
        s = zeroDeltaFun()        
    end
    
end    


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Functions implemented in this file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

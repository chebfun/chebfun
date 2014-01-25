classdef (InferiorClasses = {?bndfun, ?unbndfun}) deltafun < fun
%DELTAFUN   Class for distributions based on Dirac-delta functions on arbitrary
%   intervals.
%
%   Class for approximating generalized functions on the interval [a, b].
%   The smooth or classical part of the function is approximated by a
%   CLASSICFUN object while the Dirac delta functions are represented by a
%   pair of properties: DELTAMAG and LOCATION, which store the magnitude
%   and the location of the delta functions respectively. DELTAMAG is
%   generally a matrix, with its first row representing the delta functions,
%   while derivativeds of delta functions are represented by higher rows of
%   this matrix. LOCATION is a vector, each element of which corresponds to
%   the location of a column in the DELTAMAG matrix.
%
% Constructor inputs:
%   DELTAFUN(IMPULSES, LOCATION) creates a DELTAFUN object with an empty
%   FUNPART, while the delta functions and their locations are specified by
%   DELTAMAG and LOCATION.
% 
%   DELTAFUN(FUNPART, DELTAMAG, LOCATION) creates a DELTAFUN object with
%   FUNPART as its smooth function, while DELTAMAG and LOCATION specify the
%   delta functions in this object.
%
%   DELTAFUN(FUNPART, DELTAMAG, LOCATION, PREF) is the same as above but
%   uses PREF to pass any preferences.
%
% See also CLASSICFUN, ONEFUN, FUN

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

    %% Properties of DELTAFUN objects
    properties ( Access = public )
        % Smooth part of the representation.
        funPart     % (classical function which is a CLASSICFUN object)
        
        % DELTAMAG is a matrix storing signed magnitude of the delta functions
        % contained in a DELTAFUN object. The first row corresponds to the
        % delta functions themselves while higher order rows represent
        % derivatives of delta functions.
        deltaMag % double matrix [M X N]
        
        % Location of the delta functions.
        location % double vector [1 X N]            
    end
    
    %% DELTAFUN CLASS CONSTRUCTOR:
    methods ( Static = true )
        function obj = deltafun(funPart, deltaMag, location, pref)
            %%
            % Check for preferences in the very beginning.
            if ( (nargin < 4) || isempty(pref) )
                % Determine preferences if not given.
                pref = chebpref();
            else
                % Merge if some preferences are given.
                pref = chebpref(pref);
            end
                        
            %% Cases based on the number of arguments
            % Case 0: No input arguments, return an empty object.
            if ( nargin == 0 )
                obj.funPart = [];
                obj.deltaMag = [];
                obj.location = [];                               
                return
            end
            
            %%
            % Case 1: One input argument.
            % This input should be a fun.
            if ( nargin == 1 )
                if ( isempty(funPart) )
                    obj.funPart = [];
                elseif ( ~isa(funPart, 'fun') )
                    error( 'DELTAFUN:ctor', 'funPart must be a fun' );
                else
                    obj.funPart = funPart;
                end
                deltaMag = [];
                location = [];
            end
            %%
            % Case 2: Two input arguments.
            % Assume the argumenst passed are impulses and their locations. 
            if ( nargin == 2 )
                % Assign empty fun:
                obj.funPart = [];
                location = deltaMag;
                deltaMag = funPart;
                % Do no checks here, they are all done below.
            end
            
            %%
            % Case 3: Three input arguments.
            if ( nargin >= 3)                            
                if ( isempty(funPart) )
                    obj.funPart = [];
                elseif ( ~isa(funPart, 'fun') )
                    error( 'DELTAFUN:ctor', 'funPart must be a fun' );
                else
                    obj.funPart = funPart;
                end
                % Do no checks here, they are all done below.
            end
               
            %% Check all the arguments:            
            
            % If one of impulses or location is empty, make both empty:
            if ( isempty(deltaMag) || isempty(location) )
                deltaMag = [];
                location = [];    
            end            
            
            % Make sure location is a row vector:
            if ( ~isempty(location) )
                if ( min(size(location)) > 1 )
                    error('DELTAFUN:dim', 'location should be a vector');
                end
                location = location(:).';
            end
            
            % Check sizes:
            if ( ~isempty(deltaMag) )
                if ( size(deltaMag, 2) ~= length(location) )
                    error('DELTAFUN:dim', 'Impulse matrix should have the same number of columns as locations' );
                end
            end                         
      
            % Locations of delta functions should be within the domain:
            if ( ~isempty(location) )
                if ( ~isempty(obj.funPart) )
                    dom = obj.funPart.domain;
                    if( max(location) > dom(2) || min(location) < dom(1)  )
                        error('DELTAFUN:domain', 'Location of a delta fun is outside the domain');
                    end
                end
            end
            
            % All checks done, assign inputs to the current object:
            obj.deltaMag = deltaMag;
            obj.location = location;
                                
            % Simplify to merge redundant impulses:
            obj = simplify(obj, pref);
        end
    end
    
    
    %% METHODS IMPLEMENTED BY THIS CLASS.
    methods
        
        % True if the DELTAFUN object has no delta functions       
        out = anyDelta(f)

        % Compose a DELTAFUN with an operator or another DELTAFUN
        f = compose(f, op, g, pref)
        
        % Complex conjugate of a DELTAFUN.
        f = conj(f)
        
        % DELTAFUN obects are not transposable.
        f = ctranspose(f)
        
        % DELTAFUN display of useful data
        data = dispData(f)
        
        % Indefinite integral of a DELTAFUN.
        [F, jumpVals, locations] = cumsum(f)
        
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
        
        % Compute the inner product of two DELTAFUN objects.
        out = innerProduct(f, g)
                
        % True for an empty DELTAFUN.
        out = isempty(f)
        
        % Test if DELTAFUN objects are equal.
        out = isequal(f, g)
        
        % Test if a DELTAFUN is bounded.
        out = isfinite(f)
        
        % Test if a DELTAFUN is unbounded.
        out = isinf(f)
        
        % Test if a DELTAFUN has any NaN values.
        out = isnan(f)
        
        % True for real DELTAFUN.
        out = isreal(f)
                
        % Ture if the DELTAFUN object is zero
        out = iszero(f)
        
        % Length of a DELTAFUN.
        len = length(f)
        
        % Global maximum of a DELTAFUN on [-1,1].
        [maxVal, maxPos] = max(f)
        
        % Global minimum of a DELTAFUN on [-1,1].
        [minVal, minPos] = min(f)
        
        % Global minimum and maximum of a DELTAFUN on [-1,1].
        [vals, pos] = minandmax(f)
        
        % Subtraction of two DELTAFUN objects.
        f = minus(f, g)
        
        % Left matrix divide for DELTAFUN objects.
        X = mldivide(A, B)
        
        % Right matrix divide for a DELTAFUN.
        X = mrdivide(B, A)
        
        % Multiplication of DELTAFUN objects.
        f = mtimes(f, c)
        
        % Basic linear plot for DELTAFUN objects.
        varargout = plot(f, varargin)
        
        % Data for plotting a FUN
        data = plotData(f, g);

               
        % Addition of two DELTAFUN objects.
        f = plus(f, g)       
        
        % Dividing two DELTAFUNs
        f = rdivide(f, g)
        
        % Real part of a DELTAFUN.
        f = real(f)
        
        % Restrict a SINGFUN to a subinterval.
        f = restrict(f, s)
        
        % Roots of a DELTAFUN.
        out = roots(f, varargin)
        
        % Size of a DELTAFUN.
        [siz1, siz2] = size(f, varargin)
        
        % Simplify a DELTAFUN
        f = simplify(f, pref)
        
        % Definite integral of a DELTAFUN.
        out = sum(f, dim)
        
        % DELTAFUN multiplication.
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
        [A, v] = cleanColumns(A, v, pref);
        
        % Create map
        map = createMap(ends)
        
        % remove zero trailing rows
        A = cleanRows(A, pref);
        
        % Find intersection based on some tolerance
        [x, idxV, idxW] = numIntersect( V, W, tol)
        
        % Constructor shortcut
        f = make(varargin)
        % Merge columns of a matrix based on duplicate values in v.
        [A, v] = mergeColumns(A, v, pref)
        
        % Merge impulse matrix
        [D, w] = mergeImpulses(A, v, B, u);
        
        % Costruct a zero DELTAFUN
        s = zeroDeltaFun(domain)
        
    end
    
end
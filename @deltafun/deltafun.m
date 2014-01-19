classdef (InferiorClasses = {?bndfun, ?unbndfun}) deltafun < fun
    %DELTAFUN   Class for distributions based on Dirac-delta functions on arbitrary
    %   intervals.
    %
    %   Class for approximating generalized functions on the interval [a, b].
    %
    %   DELTAFUN class description
    %   [TODO]:
    %
    %   [TODO]: Calling Sequence
    %
    % See also PREF
    
    % Copyright 2013 by The University of Oxford and The Chebfun Developers.
    % See http://www.chebfun.org/ for Chebfun information.
    
    %% Properties of DELTAFUN objects
    properties ( Access = public )
        % Smooth part of the representation.
        funPart     % (classical function which is a CLASSICFUN object)
        
        % [TODO]: Change this documentation:
        % IMPULSES is a three-dimensional array storing information about the
        % values of the CHEBFUN object at the points in DOMAIN. The rows
        % correspond to the breakpoints in the DOMAIN vector, and if M > 1 then
        % the columns correspond to the columns in an array-valued CHEBFUN.
        % Thus, F.IMPULSES(:, :, 1) is a matrix consisting of the values of
        % each column of F at each breakpoint. The third dimension is used for
        % storing information about higher-order delta functions that may be
        % present at breakpoints. (See "help dirac" for more details.)
        impulses
        
        % location
        location               
    end
    
    %% DELTAFUN CLASS CONSTRUCTOR:
    methods ( Static = true )
        function obj = deltafun(funPart, impulses, location, pref)
            %%
            % Check for preferences in the very beginning.
            if ( (nargin < 4) || isempty(pref) )
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
                obj.impulses = [];
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
                impulses = [];
                location = [];
            end
            %%
            % Case 2: Two input arguments.
            % Assume the argumenst passed are impulses and their locations. 
            if ( nargin == 2 )
                % Assign empty fun:
                obj.funPart = [];
                location = impulses;
                impulses = funPart;
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
            if ( isempty(impulses) || isempty(location) )
                impulses = [];
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
            if ( ~isempty(impulses) )
                if ( size(impulses, 2) ~= length(location) )
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
            
            % All checks done, assign inputs to object:
            obj.impulses = impulses;
            obj.location = location;
                                
            % Simplify to merge redundant impulses:
            obj = simplify(obj);
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
        
        % Indefinite integral of a DELTAFUN.
        [F, jumpVals, locations] = cumsum(f)
        
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
        
        % Compute the inner product of two DELTAFUN objects.
        out = innerProduct(f, g)
        
        % Innerproduct, equivalent to action of a distribution 
        % on a FUN
        out = ip(f,g)
        
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
        
        % Definite integral of a DELTAFUN.
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
        
        % Create map
        map = createMap(ends)
        
        % remove zero trailing rows
        A = cleanRows(A);
        
        % Find intersection based on some tolerance
        [x, idxV, idxW] = numIntersect( V, W, tol)
        
        % Constructor shortcut
        %[TODO]: revisit this
        f = make(varargin)
        % Merge columns of a matrix based on duplicate values in v.
        [A, v] = mergeColumns(A, v)
        
        % Merge impulse matrix
        [D, w] = mergeImpulses(A, v, B, u);
        % Retrieve and modify preferences for this class.
        prefs = pref(varargin)
        
        % Costruct a zero DELTAFUN
        s = zeroDeltaFun(domain)
        
    end
    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Functions implemented in this file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

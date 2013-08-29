classdef singfun    
%SINGFUN   Class for functions with singular endpoint behavior.
%
%   Class for approximating singular functions on the interval [-1,1] using
%   a smooth part with no singularities and two singular factors (1+x).^a
%   and (1-x).^b, where a and b are negative real number or positive fractions. 
%
%   SINGFUN class description
%
%   The singfun class represents a function of the form 
%
%          f(x) = s(x) (1+x)^a (1-x)^b 
%
%   on the interval [-1,1]. The exponents a and b are assumed
%   to be real. The constructor is supplied with a handle that evaluates the 
%   function f at any given points within the interval [-1, 1]. The endpoint 
%   values are likely to return Inf or NaN results.
%
%   Ideally, the "smooth" function s is analytic, or at least much more
%   compactly represented than f is. The resulting object can be used to
%   evaluate and operate on the function f. If a and b are unknown at the 
%   time of construction, the constructor will try to determine appropriate 
%   values automatically by sampling the function handle. Note, however, that 
%   this process is not completely robust, and the singularity terms in 
%   general do not perfectly factor out singular behavior. The constructor 
%   can be forced to consider only integer exponents.
%
%   Multiplication and division are as good as the corresponding operations
%   on the smooth part. Addition and subtraction are much less reliable, as
%   the sum of two singfuns with different exponents is not necessarily a
%   singfun, nor a smooth function. If all but integer exponents can be
%   factored out of the summands, the process is fine, but in other
%   circumstances the process may throw an error.
%
% See also SINGFUN.PREF

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

    %% Properties of SINFGUN objects
    properties ( Access = public )
        % Smooth part of the representation.
        smoothPart  % (smoothfun)
        
        % Exponents of the singularities at the two endpoints.
        exponents   % (1x2 double)
        
        % A cell array containing the types of singularities at the endpoints.
        singType = {};   % (1x2 cell)        
     end
    
    %% CLASS CONSTRUCTOR:
    methods ( Static = true )
        function obj = singfun(op, exponents, singType, pref) 
            %%
            % Check for preferences in the very beginning.
            % Determine preferences if not given, merge if some are given:
            if ( nargin < 4 || isempty(pref) )
                pref = singfun.pref;
            else        
                pref = singfun.pref(pref);
            end
            %%
            % Check for cases based on the number of arguments            
            
            %%
            % no input arguments: return an empty object               
            if ( nargin == 0 )   
                obj.smoothPart = [];
                obj.exponents = [];
                obj.singType = {};
                return
            end
            
            %%
            if ( nargin == 1 )
                % only operator passed, assume a fractional pole at each 
                % end point               
                exponents = [];
                obj.singType = {'sing', 'sing'};                
            end
            %%
            if ( nargin == 2 || ~isempty(exponents) )
                % exponents passed, discard the values
                % given in singType and use the
                % information given in exponents.
                obj.exponents = exponents;
                if ( nargin == 2 || isempty(singType) )
                    % if the user has not provided the type of
                    % singularity, figure it out
                    obj = classifyExponents(obj);
                else
                    % Singularity types given, make sure the strings are OK.
                    if ( ~isa(singType, 'cell') )
                        error( 'CHEBFUN:SINGFUN:constructor', ...
                               'singType must be a 1x2 cell with two strings');
                    end
                    obj.singType = singType;
                    checkSingTypes(obj);                    
                end
                              
            end
                
            %%
            if ( nargin >= 3 && isempty(exponents) )                
                if ( isempty(singType) )
                    % Singulrity types and exponents not given. Assume
                    % fractional poles or generic singularities if not given
                    obj.singType = {'sing', 'sing'};
                else
                    % Singularity types given, make sure the strings are OK.
                    if ( ~isa(singType, 'cell') )
                        error( 'CHEBFUN:SINGFUN:constructor', ...
                               'singType must be a 1x2 cell with two strings');
                    end
                    obj.singType = singType;
                    checkSingTypes(obj);
                end
            end                        
            
            %%
            % Determine and factor out singular terms if exponents 
            % are not given
            if ( isempty(obj.exponents) )
                obj.exponents = singfun.findSingExponents(op, obj.singType);
                % update SINGTYPE based on EXPONENTS
                obj = classifyExponents(obj);               
            end
               
            % update the operator based on the values in exponents.
            smoothOp = singfun.singOp2SmoothOp(op, obj.exponents);                        
            obj.smoothPart = singfun.constructSmoothPart(smoothOp, pref);
        end
    end
    
    %% 

    
    %% METHODS IMPLEMENTED BY THIS CLASS.
    methods
        
        % Complex conjugate of a SINGFUN.
        f = conj(f)                
        
        % Classify epxonents of a SINGFUN
        f = calssifyExponents(f)
        
        % Check the strings classifying singTypes
        out = checkSingTypes(f)
        
        % SINGFUN obects are not transposable.
        f = ctranspose(f)

        % Indefinite integral of a SINGFUN.
        f = cumsum(f, m, pref)

        % Derivative of a SINGFUN.
        f = diff(f, k)
        
        % Evaluate a SINGFUN.
        y = feval(f, x)

        % Flip columns of an array-valued SINGFUN object.
        f = fliplr(f)
        
        % Flip/reverse a SINGFUN object.
        f = flipud(f)
        
        % Imaginary part of a SINGFUN.
        f = imag(f)
     
        % True for an empty SINGFUN.
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
        s = constructSmoothPart( op, pref)
        
        % method for finding the order of singularities
        exponents = findSingExponents( op, singType )
        
        % method for finding integer order singularities, i.e. poles
        poleOrder = findPoleOrder( op, SingEnd)
        
        % method for finding fractional order singularities (+ve or -ve).
        barnchOrder = findSingOrder( op, SingEnd)
        
        % method for converting a singular op to a smooth op
        op = singOp2SmoothOp( op, exponents, tol )
        
        % Retrieve and modify preferences for this class.
        prefs = pref(varargin)
        
        % Costruct a zero SINGFUN
        s = zeroSingFun()        
    end
    
end    
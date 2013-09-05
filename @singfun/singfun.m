classdef singfun    
%SINGFUN   Class for functions with singular endpoint behavior.
%
%   Class for approximating singular functions on the interval [-1,1] using a
%   smooth part with no singularities and two singular factors (1+x).^a and
%   (1-x).^b, where a and b are real numbers.
%
%   SINGFUN class description
%
%   The SINGFUN class represents a function of the form 
%
%          f(x) = s(x) (1+x)^a (1-x)^b 
%
%   on the interval [-1,1]. The exponents a and b are assumed to be real. The
%   constructor is supplied with a handle that evaluates the function f at any
%   given points within the interval [-1,1]. The endpoint values are likely to
%   return Inf or NaN results.
%
%   Ideally, the "smooth" function s is analytic, or at least much more
%   compactly represented than f is. The resulting object can be used to
%   evaluate and operate on the function f. If a and b are unknown at the time
%   of construction, the constructor will try to determine appropriate values
%   automatically by sampling the function handle. Note, however, that this
%   process is not completely robust, and the singular terms in general do not
%   perfectly factor out singular behavior. The constructor can be forced to
%   consider only integer exponents.
%
%   [TODO]: Describe calling sequence.
%
%   Multiplication and division are as good as the corresponding operations on
%   the smooth part. Addition and subtraction are much less reliable, as the sum
%   of two SINGFUN objects with different exponents is not necessarily a
%   SINGFUN, nor a smooth function. If all but integer exponents can be factored
%   out of the summands, the process is fine, but in other circumstances the
%   process may throw an error.
%
% See also PREF

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

    %% Properties of SINGFUN objects
    properties ( Access = public )
        % Smooth part of the representation.
        smoothPart      % (smoothfun)
        
        % Exponents of the singularities at the two endpoints.
        exponents       % (1x2 double)
                
        % [TODO]: Shold exponentTol be a property?
     end
    
    %% CLASS CONSTRUCTOR:
    methods ( Static = true )
        function obj = singfun(op, exponents, singType, pref) 
            %%
            % Check for preferences in the very beginning.            
            if ( (nargin < 4) || isempty(pref) )
                % Determine preferences if not given.
                pref = singfun.pref;
            else
                % Merge if some preferences are given.
                pref = singfun.pref(pref);
            end
            
            
            %% Cases based on the number of arguments                        
            % Case 0: No input arguments, return an empty object.
            if ( nargin == 0 )   
                obj.smoothPart = [];
                obj.exponents = [];
                return
            end
            
            %%
            % Case 1: One input argument.
            if ( nargin == 1 )
                % Make sure the exponents are empty.
                exponents = [];                
            end
            %%
            if ( (nargin == 2) || ~isempty(exponents) )
                % Exponents passed, store them.
                obj.exponents = exponents;                                           
            end
                
            %%
            if ( nargin >= 3 && isempty(exponents) )                
                if ( isempty(singType) )
                    % Singulrity types and exponents not given. Assume
                    % fractional poles or generic singularities if not given
                    singType = {'sing', 'sing'};
                else
                    % Singularity types given, make sure the strings are OK.
                    if ( ~isa(singType, 'cell') )
                        error( 'CHEBFUN:SINGFUN:constructor', ...
                               'singType must be a 1x2 cell with two strings');
                    end                   
                    checkSingTypes(singType);
                end
            end    
            
            %%
            % Make sure that op is a funciton handle
            if ( ~isa(op, 'function_handle') )
                error( 'CHEBFUN:SINGFUN:constructor', 'First argument must be a function handle.');
            end
            
            %%
            % Determine and factor out singular terms if exponents 
            % are not given
            if ( isempty(obj.exponents) )
                obj.exponents = singfun.findSingExponents(op, singType);                
            end
               
            % Factor out singular terms from the operator based on the values
            % in exponents.
            smoothOp = singOp2SmoothOp(op, obj.exponents);
            % Construct the smooth part.
            obj.smoothPart = singfun.constructSmoothPart(smoothOp, pref);
        end
    end
    
    %% 

    
    %% METHODS IMPLEMENTED BY THIS CLASS.
    methods
        
        % Complex conjugate of a SINGFUN.
        f = conj(f)                
               
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
        
        % Retrieve and modify preferences for this class.
        prefs = pref(varargin)
        
        % Costruct a zero SINGFUN
        s = zeroSingFun()        
    end
    
end    


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Functions implemented in this file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function out = checkSingTypes(singType)
%CHECKSINGTYPES   Function to check types of exponents in a SINGFUN object.
%   The valid types can be 'sing', 'pole', 'root' or 'none'. If the type is
%   different than these four strings (ignoring case), an error message is
%   thrown.
%

%%
out(1) = any(strcmpi(singType{1}, {'pole', 'sing', 'root', 'none'}));
out(2) = any(strcmpi(singType{2}, {'pole', 'sing', 'root', 'none'}));

if ( ~all(out) )
    error('CHEBFUN:SINGFUN:checkSingTypes', 'Unknown singularity type.');
end

end

function op = singOp2SmoothOp(op, exponents)
%SINGOP2SMOOTHOP   Converts a singular operator to a smooth one by removing the 
%   singularity(ies). 
%   SINGOP2SMOOTHOP(OP, EXPONENTS) returns a smooth operator OP by removing the
%   singularity(ies) at the endpoints -1 and 1. EXPONENTS are the order of the
%   singularities. 
%
%   For examples, opNew = singOp2SmoothOp(opOld, [-a -b]) means 
%   opNew = opOld.*(1+x).^a.*(1-x).^b
%
% See also SINGFUN.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( all(exponents) )
    % Both exponents are non trivial:
    op = @(x) op(x)./((1 + x).^(exponents(1)).*(1 - x).^(exponents(2)));
    
elseif ( exponents(1) )
    % (1+x) factor at the left end point:
    op = @(x) op(x)./(1 + x).^(exponents(1));
    
elseif ( exponents(2) )
    % (1-x) factor at the right end point:
    op = @(x) op(x)./(1 - x).^(exponents(2));
    
end

end




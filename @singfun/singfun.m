classdef (InferiorClasses = {?chebtech2, ?chebtech1}) singfun < onefun %(See Notes)
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
%   Ideally, the "smooth" function s(x) is analytic, or at least much more
%   compactly represented than f is. The resulting object can be used to
%   evaluate and operate on the function f. If a and b are unknown at the time
%   of construction, the constructor will try to determine appropriate values
%   automatically by sampling the function handle. Note, however, that this
%   process is not completely robust, and the singular terms in general do not
%   perfectly factor out singular behavior. The constructor can be forced to
%   consider only integer exponents.
%
%   Multiplication and division are as good as the corresponding operations on
%   the smooth part. Addition and subtraction are much less reliable, as the sum
%   of two SINGFUN objects with different exponents is not necessarily a
%   SINGFUN, nor a smooth function. If all but integer exponents can be factored
%   out of the summands, the process is fine, but in other circumstances the
%   process may throw an error.
%
%   SINGFUN(OP) constructs a SINGFUN object from the function handle OP. It
%   first tries to determine the type and the order of the singularities or the
%   roots at the endpoints by sampling near -1 and 1 and then forms a new
%   operator which governs the smooth part of OP by factoring out the end point
%   singularities. Finally it constructs an approximation of the smooth part by
%   calling the SMOOTHFUN constructor. The type and the order of the
%   singularities along with the SMOOTHFUN representation of the smooth part are
%   stored in corresponding member fields of the instantiation.
%
%   SINGFUN(OP, [], SINGTYPE) constructs a SINGFUN object as above. SINGTYPE 
%   is 1x2 cell array of strings and type of the singularities specified by 
%   these strings may help the singularity detector to determine the order 
%   of the singularities more efficiently and save some computing time. Valid
%   strings for SINGTYPE are 'none', 'pole', 'sing' or 'root'. Note that a 
%   place holder must be given next to OP for the constructor to work 
%   properly.
%
%   SINGFUN(OP, EXPONENTS) constructs a SINGFUN object. Instead of determining
%   the singularity types and orders by sampling the function values at -1 and
%   1, the constructor takes the values saved in the 1x2 vector EXPONENTS as the
%   orders of the singularities.
%
%   SINGFUN(OP, EXPONENTS, SINGTYPE) and SINGFUN(OP, EXPONENTS, {}) are
%   equivalent to SINFGUN(OP, EXPONENTS).
%
%   SINGFUN(OP, EXPONENTS, SINGTYPE, VSCALE, HSCALE) constructs a SINGFUN 
%   object. When the smooth part s(x) is approximated, the vertical scale VSCALE
%   and the horizontal scale HSCALE are passed to the SMOOTHFUN constructor to
%   facilitate the construction. Note that any of or both of VSCALE and HSCALE 
%   can be omitted or empty.
%
%   SINGFUN(OP, EXPONENTS, SINGTYPE, VSCALE, HSCALE, PREF) constructs a SINGFUN 
%   using the preferences given by PREF. Note that any of or all of EXPONENTS, 
%   SINGTYPE, VSCALE, and HSCALE can be omitted or empty in this calling
%   sequence in the presence of PREF.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

%% NOTES:
%   It is a Matlab requirement to specify exactly which classes are inferior to
%   a given class. One might think that writing "InferiorClasses = {?smoothfun}"
%   should be OK but it turns out that subclasses do not inherit the attribute
%   of being inferior and all inferior clases/subclasses should be mentioned 
%   explicitly.

    %% Properties of SINGFUN objects
    properties ( Access = public )
        % Smooth part of the representation.
        smoothPart      % (SMOOTHFUN)
        
        % Exponents of the singularities at the two endpoints.
        exponents       % (1x2 double)
    end
    
    %% CLASS CONSTRUCTOR (IMPLEMENTED BY THIS M-FILE):
    methods
        function obj = singfun(op, exponents, singType, vscale, hscale, pref)   
            
            %% Get Preferences
            % Check for preferences in the very beginning.
            if ( (nargin < 6) || isempty(pref) )
                % Determine preferences if not given.
                pref = chebfunpref();
            else
                % Merge if some preferences are given.
                pref = chebfunpref(pref);
            end           
            
            %% Cases based on the number of arguments
            % Case 0: No input arguments, return an empty object.
            if ( nargin == 0 )
                obj.smoothPart = [];
                obj.exponents = [];
                return
            end
            
            % Case 1: One input argument.
            if ( nargin == 1 )
                if ( isa(op, 'singfun') )
                    obj = op;
                    return
                elseif ( isa(op, 'smoothfun') )
                    % if OP is a SMOOTHFUN, cast it to a SINGFUN:
                    obj.smoothPart = op;
                    obj.exponents = [0, 0];
                    
                    return
                else
                    % Make sure the exponents are empty.
                    exponents = [];
                end
            end
            
            % Case 2: Two input arguments.
            if ( (nargin == 2) || ~isempty(exponents) )
                % Exponents passed, store them.
                obj.exponents = exponents;
            end
                        
            % Case 3: Three or more input arguments.
            % The user can choose a singularity detection algorithm by passing
            % appropriate strings in the argument "singType", which is a cell
            % array. If, however, the user doesn't provide any preferences
            % regarding the algorithm, the most generic algorithm, which tries
            % to find fractional order singularities is used.
            if ( nargin >= 3 && isempty(exponents) )
                if ( isempty(singType) )
                    % Singularity types and exponents not given. Assume
                    % fractional poles or generic singularities at each
                    % end.
                    singType = {'sing', 'sing'};
                else
                    % Singularity types given, make sure the strings are OK.
                    checkSingTypes(singType);
                end
            else
                if ( isempty(exponents) )
                    singType = {'sing', 'sing'};
                end
            end
            
            %% Misc Checks on the Inputs
            
            % Make sure that op is a function handle or a smoothfun:
            if ( ~isa(op, 'function_handle') && ~isa(op, 'smoothfun') )
                error( 'CHEBFUN:SINGFUN:constructor', ...
                    'First argument must be a function handle or a smoothfun.');
            end
            
            % Check to avoid array-valued operators:
            if ( size(feval(op, 0), 2) > 1 )
                error( 'CHEBFUN:SINGFUN:constructor', ...
                    'SINGFUN class does not support array-valued objects.' );
            end
            
            % Handle cases for vscale and hscale.
            if ( nargin <= 3 )
                % If both are not passed as argument, set them to empty.
                vscale = [];
                hscale = [];
            end
            
            if ( nargin == 4 )
                % vscale passed, but hscale missing, set it to empty:
                hscale = [];
            end                   
            
            %% Find Exponents
            % If exponents were passed, make sure they are in correct shape.
            if ( ~isempty(exponents) )
                if ( any(size(exponents) ~= [1, 2]) || ...
                        ~isa(exponents, 'double') )
                    error( 'CHEBFUN:SINGFUN:constructor', ...
                        'Exponents must be a 1X2 vector of doubles.' );
                end
            else
                % Exponents not given, determine them.
                obj.exponents = singfun.findSingExponents(op, singType);
            end
            
            % If a smoothfun has been passed as the op, store it directly:
            if ( isa(op, 'smoothfun') )
                obj.smoothPart = op;
            else
                %% Construct New Function Handle
                % Factor out singular terms from the operator based on the values
                % in exponents.
                smoothOp = singOp2SmoothOp(op, obj.exponents);
                
                % Construct the smooth part.
                obj.smoothPart = singfun.constructSmoothPart(smoothOp, ...
                    vscale, hscale, pref);
            end
        end
    end
    
    %% METHODS (NON-STATIC) IMPLEMENTED BY THIS CLASS.
    methods
        
        % SINGFUN logical AND.
        h = and(f, g)

        % True if any element of a SINGFUN is a nonzero number, ignoring NaN.
        a = any(f, dim)

        % Cancel the negative exponents of a SINGFUN.
        f = cancelExponents(f)

        % Convert an array of ONEFUN objects into an array-valued ONEFUN.
        f = cell2mat(f)
        
        % Complex conjugate of a SINGFUN.
        f = conj(f)
        
        % SINGFUN objects are not transposable.
        f = ctranspose(f)
        
        % Indefinite integral of a SINGFUN.
        f = cumsum(f, dim)
        
        % Derivative of a SINGFUN.
        f = diff(f, k, dim)
        
        % Extract information for DISPLAY.
        info = dispData(f)
        
        % Extract roots at the boundary points -1 and 1.
        [f, rootsLeft, rootsRight] = extractBoundaryRoots(f, numRoots)
        
        % Evaluate a SINGFUN.
        y = feval(f, x)
        
        % SINGFUN does not support FIX.
        g = fix(f);
        
        % Flip columns of an array-valued SINGFUN object.
        f = fliplr(f)
        
        % Flip/reverse a SINGFUN object.
        f = flipud(f)
        
        % SINGFUN does not support FLOOR.
        g = floor(f);
        
        % Get method:
        val = get(f, prop);
        
        % Imaginary part of a SINGFUN.
        f = imag(f)
        
        % Compute the inner product of two SINGFUN objects.
        out = innerProduct(f, g)
        
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
        
        % True for smooth SINGFUN.
        out = issmooth(f)
        
        % True for zero SINGFUN objects
        out = iszero(f)
        
        % Length of a SINGFUN.
        len = length(f)
        
        % SINGFUN logical.
        f = logical(f)
        
        % Convert an array-valued SINGFUN into an ARRAY of SINGFUN objects.
        g = mat2cell(f, M, N)
        
        % Global maximum of a SINGFUN on [-1,1].
        [maxVal, maxPos] = max(f)
        
        % Global minimum of a SINGFUN on [-1,1].
        [minVal, minPos] = min(f)
        
        % Global minimum and maximum of a SINGFUN on [-1,1].
        [vals, pos] = minandmax(f)
        
        % Subtraction of two SINGFUN objects.
        f = minus(f, g)
        
        % Left matrix divide for ONEFUN objects.
        X = mldivide(A, B)

        % Right matrix divide for a ONEFUN.
        X = mrdivide(B, A)
        
        % Multiplication of SINGFUN objects.
        f = mtimes(f, c)
        
        % Compute a Legendre series expansion of a ONEFUN object:
        c = legpoly(f)

        % Estimate of the norm of a SINGFUN object.
        out = normest(f)
        
        % SINGFUN logical NOT.
        f = not(f)

        % SINGFUN logical OR.
        h = or(f, g)
        
        % Basic linear plot for SINGFUN objects.
        varargout = plot(f, varargin)
        
        % Plot3() for SINGFUN objects.
        varargout = plot3(f, g, h)
        
        % Obtain data used for plotting a SINGFUN object.
        data = plotData(f, g, h)
        
        % Addition of two SINGFUN objects.
        f = plus(f, g)
        
        % Polynomial coefficients of a ONEFUN.
        out = poly(f)
        
        % SINGFUN power function.
        f = power(f, b)

        % QR factorisation of an array-valued ONEFUN.
        [f, R, E] = qr(f, flag, methodFlag)
        
        % Dividing two SINGFUNs.
        f = rdivide(f, g)
        
        % Real part of a SINGFUN.
        f = real(f)
        
        % Simplify the exponents of a SINGFUN.
        f = simplifyExponents(f)

        % Restrict a SINGFUN to a subinterval.
        f = restrict(f, s)
        
        % Roots of a SINGFUN in the interval [-1,1].
        out = roots(f, varargin)
        
        % Convert a SINGFUN to a SMOOTHFUN.
        f = singFun2SmoothFun(f) 
        
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
        % SmoothPart constructor
        s = constructSmoothPart( op, vscale, hscale, pref )
        
        % Method for finding the order of singularities
        exponents = findSingExponents( op, singType )
        
        % Find integer order singularities, i.e. poles
        poleOrder = findPoleOrder( op, SingEnd )
        
        % Finding fractional order singularities (+ve or -ve).
        branchOrder = findSingOrder( op, SingEnd )
        
        % Make a SINGFUN (constructor shortcut):
        f = make(varargin);
        
        % Retrieve and modify preferences for this class.
        prefs = pref(varargin)
      
        % Convert a SMOOTHFUN to a SINGFUN.
        f = smoothFun2SingFun(f) 
        
        % Construct a zero SINGFUN
        s = zeroSingFun()
    end
    
end

%% OTHER FUNCTIONS IMPLEMENTED IN THIS M-FILE:

function out = checkSingTypes(singType)
%CHECKSINGTYPES   Function to check types of exponents in a SINGFUN object.
%   The valid types can be 'sing', 'pole', 'root' or 'none'. If the type is
%   different than these four strings (ignoring case), an error message is
%   thrown.

if ( ~isa(singType, 'cell') )
    error( 'CHEBFUN:SINGFUN:constructor', ...
        'singType must be a 1x2 cell with two strings');
end

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
%   For example, opNew = singOp2SmoothOp(opOld, [-a -b]) means
%   opNew = opOld.*(1+x).^a.*(1-x).^b
%
% See also SINGFUN.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
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

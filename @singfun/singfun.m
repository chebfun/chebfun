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
%   SINGFUN(OP, DATA) constructs a SINGFUN object using the data in the MATLAB
%   structure DATA.  DATA fields recognized by SINGFUN are.
%     DATA.SINGTYPE    (Default:  Determined by CHEBFUNPREF)
%         A 1x2 cell array of strings and type of the singularities specified
%         by these strings may help the singularity detector to determine the
%         order of the singularities more efficiently and save some computing
%         time.  Valid strings for SINGTYPE are 'none', 'pole', 'sing' or
%         'root'.
%     DATA.EXPONENTS   (Default:  Empty)
%         If DATA.EXPONENTS is nonempty, then instead of determining the
%         singularity types and orders by sampling the function values at -1
%         and 1, the constructor takes the values saved in the 1x2 vector
%         EXPONENTS as the orders of the singularities.
%   If any fields in DATA are empty or not supplied, or if DATA itself is empty
%   or not supplied, appropriate default values are set.  Any fields in DATA
%   which are not recognized will be passed as-is to the SMOOTHFUN constructor.
%
%   SINGFUN(OP, DATA, PREF) constructs a SINGFUN using the preferences given by
%   PREF.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DEVELOPER NOTES:
%   It is a Matlab requirement to specify exactly which classes are inferior to
%   a given class. One might think that writing "InferiorClasses = {?smoothfun}"
%   should be OK but it turns out that subclasses do not inherit the attribute
%   of being inferior and all inferior clases/subclasses should be mentioned 
%   explicitly.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CLASS PROPERTIES:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    properties ( Access = public )
        % Smooth part of the representation.
        smoothPart      % (SMOOTHFUN)
        
        % Exponents of the singularities at the two endpoints.
        exponents       % (1x2 double)
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CLASS CONSTRUCTOR:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods ( Access = public, Static = false )
        
        function obj = singfun(op, data, pref)
            
            % Parse inputs.
            if ( nargin < 1 )
                % No input arguments; return an empty SINGFUN.
                obj.smoothPart = [];
                obj.exponents = [];
                return
            end

            if ( (nargin < 2) || isempty(data) )
                    data = struct();
            end

            if ( (nargin < 3) || isempty(pref) )
                pref = chebfunpref();
            else
                pref = chebfunpref(pref);
            end

            data = parseDataInputs(data, pref);

            % Leave SINGFUN OPs alone, and upgrade SMOOTHFUN OPs to SINGFUNs.
            if ( nargin == 1 )
                if ( isa(op, 'singfun') )
                    obj = op;
                    return
                elseif ( isa(op, 'smoothfun') )
                    obj.smoothPart = op;
                    obj.exponents = [0, 0];
                    return
                end
            end

            % Find exponents.poly
            if ( ~isempty(data.exponents) )
                % Check values of supplied exponents.
                if ( any(size(data.exponents) ~= [1, 2]) || ...
                        ~isa(data.exponents, 'double') )
                    error('CHEBFUN:SINGFUN:singfun:badExponents', ...
                        'Exponents must be a 1x2 vector of doubles.');
                end

                % If any of EXPONENTS is NaN, try to determine:
                maskNaN = isnan(data.exponents);
                if ( any(maskNaN) )
                    tmpExps = singfun.findSingExponents(op, data.singType);
                    data.exponents(maskNaN) = tmpExps(maskNaN);
                end

                obj.exponents = data.exponents;
            else
                % Exponents not given.  Determine them.
                obj.exponents = singfun.findSingExponents(op, data.singType);
            end

            % Make sure that op is a function handle or a smoothfun.
            if ( ~isa(op, 'function_handle') && ~isa(op, 'smoothfun') )
                error( 'CHEBFUN:SINGFUN:singfun:badOp', ...
                    ['First argument must be a function handle or a ', ...
                     'SMOOTHFUN, not a %s.'], class(op));
            end

            % Check to avoid array-vapolylued operators.
            if ( size(feval(op, 0), 2) > 1 )
                error('CHEBFUN:SINGFUN:singfun:arrayValued', ...
                    'SINGFUN does not support array-valued construction.');
            end

            % Extrapolate when the given function blows up.
            if ( any(obj.exponents < 0) )
                pref.techPrefs.extrapolate = true;
            end

            % Construct the smoothPart.
            if ( isa(op, 'smoothfun') )
                % smoothPart was handed to us.
                obj.smoothPart = op;
            else
                % Construct New Function Handle
                
                % Loosen tolerance:
                if ( any(obj.exponents) )
                    pref.eps = max(pref.eps, 1e-14);
                end
                
                % Factor out singular terms from the operator based on the values
                % in exponents.
                smoothOp = singOp2SmoothOp(op, obj.exponents);
                obj.smoothPart = singfun.constructSmoothPart(smoothOp, data, ...
                    pref);
            end
        end
        
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CLASS METHODS:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods ( Access = public, Static = false )
        
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
        
        % Extract columns (or rows) of SINGFUN.
        f = extractColumns(f, colIdx)
        
        % Evaluate a SINGFUN.
        y = feval(f, x)
        
        % SINGFUN does not support FIX.
        g = fix(f)
        
        % Flip columns of an array-valued SINGFUN object.
        f = fliplr(f)
        
        % Flip/reverse a SINGFUN object.
        f = flipud(f)
        
        % SINGFUN does not support FLOOR.
        g = floor(f)
        
        % Get method:
        val = get(f, prop)
        
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
        
        function out = isPeriodicTech(f)
        %ISPERIODICTECH   Test if the smooth part of f is is constructed with a 
        %basis of periodic functions. 
        
            % Calls ISPERIODICTECH on the SMOOTHFUN part.
            out = isPeriodicTech(f.smoothPart);
        end
        
        % True for real SINGFUN.
        out = isreal(f)
        
        % True for smooth SINGFUN.
        out = issmooth(f)
        
        % True for zero SINGFUN objects
        out = iszero(f)
        
        % Legendre coeffs of a SINGFUN.
        c = legcoeffs(f, n)
        
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
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% STATIC METHODS:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods ( Access = public, Static = true )
        
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
      
        % Convert a SMOOTHFUN to a SINGFUN.
        f = smoothFun2SingFun(f) 
        
        % Construct a zero SINGFUN
        s = zeroSingFun()
        
    end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% OTHER FUNCTIONS IMPLEMENTED IN THIS M-FILE:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function checkSingTypes(singType)
%CHECKSINGTYPES   Function to check types of exponents in a SINGFUN object.
%   The valid types can be 'sing', 'pole', 'root' or 'none'. If the type is
%   different than these four strings (ignoring case), an error message is
%   thrown.

if ( ~isa(singType, 'cell') )
    error( 'CHEBFUN:SINGFUN:checkSingTypes:notCell', ...
        'singType must be a 1x2 cell with two strings');
end

check(1) = any(strcmpi(singType{1}, {'pole', 'sing', 'root', 'none'}));
check(2) = any(strcmpi(singType{2}, {'pole', 'sing', 'root', 'none'}));

if ( ~all(check) )
    error('CHEBFUN:SINGFUN:checkSingTypes:badType', ...
        'Unknown singularity type.');
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

function data = parseDataInputs(data, pref)
%PARSEDATAINPUTS   Parse inputs from the DATA structure and assign defaults.

if ( ~isfield(data, 'exponents') || isempty(data.exponents) )
    data.exponents = [];
end

if ( ~isfield(data, 'singType') || isempty(data.singType) )
    defaultSingType = pref.blowupPrefs.defaultSingType;
    data.singType = {defaultSingType, defaultSingType};
end

end

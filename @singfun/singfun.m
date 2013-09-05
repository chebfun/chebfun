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
        
        % A cell array containing the types of singularities at the endpoints.
        singType = {};  % (1x2 cell)        
        % [TODO]: What values may this take? What do they mean? Why is it needed?
        
        % [TODO]: Shold exponentTol be a property?
     end
    
    %% CLASS CONSTRUCTOR:
    methods ( Static = true )
        function obj = singfun(op, exponents, singType, pref) 
            %%
            % Check for preferences in the very beginning.
            % Determine preferences if not given, merge if some are given:
            if ( (nargin < 4) || isempty(pref) )
                pref = singfun.pref;
            else        
                pref = singfun.pref(pref);
            end
            %%
            % Check for cases based on the number of arguments            
            
            %%
            % No input arguments: return an empty object               
            if ( nargin == 0 )   
                obj.smoothPart = [];
                obj.exponents = [];
                obj.singType = {};
                return
            end
            
            %%
            if ( nargin == 1 )
                % Only operator passed, assume a fractional pole at each 
                % end point               
                exponents = [];
                obj.singType = {'sing', 'sing'};                
            end
            %%
            if ( (nargin == 2) || ~isempty(exponents) )
                % Exponents passed, discard the values
                % given in singType and use the
                % information given in exponents.
                obj.exponents = exponents;
                if ( (nargin == 2) || isempty(singType) )
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
            % Make sure that op is a funciton handle
            if ( ~isa(op, 'function_handle') )
                error( 'CHEBFUN:SINGFUN:constructor', 'First argument must be a function handle.');
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
            smoothOp = singOp2SmoothOp(op, obj.exponents);
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

function f = classifyExponents(f)
%CLASSIFYEXPONENTS   Function to assign types of exponents in a SINGFUN object.
%   Based on the values in F.EXPONENTS, this functions decides the type that 
%   should be assigned to F.SINGTYPE. The valid types can be 'sing', 'pole', 
%   'root' or 'none'. F.SINGTYPE is a 1X2 cell array and the pair of string
%   contained in this field describes the types of singularities at each end,
%   -1 or 1 of the SINGFUN F. These types have the following meaning:
%    
%      'pole' - A pole, i.e. a negative integer exponent at the 
%               corresponding end.
%      'sing' - A negative real exponent at the corresponding end. 
%               Can be an integer as well.
%      'root' - A root of fractional order at the corresponding end point.
%      'none' - No singularity at the end point.

%%
% Get the SINGFUN tolerance
% [TODO]: This should depend on scales, but what are the scales?
%         This TODO will be setteled once scales are finalised.
tol = singfun.pref.singfun.eps;

%%
% Store the exponents in the variable exps (for brevity):
exps = f.exponents;
% Loop on the left and right end point of the domain
for k = 1:2
    % If positive exponent
    if ( exps(k) >= 0 )
        if ( abs(exps(k) - round(exps(k))) < 100*tol )
            % Positive integer exponent, i.e. no singularity
            f.singType{k} = 'none';
        else
            % The function is bounded but there is a root of fractional order.
            f.singType{k} = 'root';
        end
    else        
        % Negative exponents        
        if ( exps(k) > -100*tol )
            % The exponent is negative but almost zero,
            % remove the singularity
            f.singType{k} = 'none';
        else
            % Non-trivial negative exponent
            if ( abs(exps(k) - round(exps(k))) < 100*tol )
                % Pole if integer valued exponent
                f.singType{k} = 'pole';
            else
                % A fractional pole, which we call 'sing'
                f.singType{k} = 'sing';
            end
        end
    end
end

end

function out = checkSingTypes(f)
%CHECKSINGTYPES   Function to check types of exponents in a SINGFUN object.
%   The valid types can be 'sing', 'pole', 'root' or 'none'. If the type is
%   different than these four strings (ignoring case), an error message is
%   thrown.
%

%%
out(1) = any(strcmpi(f.singType{1}, {'pole', 'sing', 'root', 'none'}));
out(2) = any(strcmpi(f.singType{2}, {'pole', 'sing', 'root', 'none'}));

if ( ~all(out) )
    error('CHEBFUN:SINGFUN:checkSingTypes', 'Unknown singularity type');
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




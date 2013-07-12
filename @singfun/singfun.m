classdef singfun
    
%SINGFUN Class for functions with singular endpoint behavior.
%   TODO: User documentation

%% SINGFUN class description
%
% The singfun class represents a function of the form 
%
%     f(x) = s(x) (1+x)^a (1-x)^b 
%
% on the interval [-1,1]. The exponents a and b are assumed
% to be real and negative. The constructor is supplied with a handle that evaluates the
% function f at any given points. However, endpoint values will not be
% sampled, due to the likelihood of Inf and NaN results.
%
% Ideally, the "smooth" function s is analytic, or at least much more
% compactly represented than f is. The resulting object can be used to
% evaluate and operate on the function f. If a and b are
% unknown at the time of construction, the constructor will try to
% determine appropriate (nonpositive) values automatically by sampling
% the function handle. Note, however, that this process is not 
% completely robust, and
% the singularity terms in general do not perfectly factor out singular
% behavior. The constructor can be forced to consider only integer
% exponents.
%
% Multiplication and division are as good as the corresponding operations
% on the smooth part. Addition and subtraction are much less reliable, as
% the sum of two singfuns with different exponents is not necessarily a
% singfun, nor a smooth function. If all but integer exponents can be
% factored out of the summands, the process is fine, but in other
% circumstances the process may throw an error.
 
% Copyright 2013 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

    %% Properties of SINFGUN objects
    properties ( Access = public )
        % Smooth part of the representation.
        smoothPart  % (smoothfun)
        
        % Exponents of the singularities at the two endpoints.
        exponents   % (1x2 double)
        
        % A cell array telling the type of singularity at the endpoints.
        singType   % (1x2 cell)
        
        % A logical array indicating which ends are singular.
        isSingEnd   % (1x2 logical)        
    end
    
    %% CLASS CONSTRUCTOR:
    methods ( Static = true )
        function obj = singfun(op, exponents, isSingEnd, singType, pref)           
            %%
            % Check for preferences in the very beginning.
            % Determine preferences if not given, merge if some are given:
            if ( nargin < 5 || isempty(pref) )
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
                obj.isSingEnd = [];
                return
            end
            
            %%
            if ( nargin == 1 )
                % only operator passed, assume a pole at each end point
                obj.isSingEnd = [1, 1];
                obj.singType = {'pole', 'pole'};                
            end
            %%
            if ( nargin == 2 || ~isempty(exponents) )
                % exponents passed, discard the values
                % given in isSingEnd, singType and use 
                % the information given in exponents.
                obj.exponents = exponents;
                tol = pref.singfun.eps;
                if ( exponents(1) < 0 && abs(exponents(1)) > 100*tol )
                    obj.isSingEnd(1) = 1;
                    if( abs(exponents(1)-round(exponents(1))) < 100*tol )
                        obj.singType{1} = 'pole';
                    else
                        obj.singType{1} = 'branch';
                    end
                else
                    obj.isSingEnd(1) = 0;
                    obj.singType{1} = 'none';
                end
                
                if ( exponents(2) < 0 && abs(exponents(2)) > 100*tol )
                    obj.isSingEnd(2) = 1;
                    if( abs(exponents(2)-round(exponents(2))) < 100*tol )
                        obj.singType{2} = 'pole';
                    else
                        obj.singType{2} = 'branch';
                    end
                else
                    obj.isSingEnd(1) = 0;
                    obj.singType{1} = 'none';                   
                end                
            end
                
            if ( nargin == 3 && isempty(exponents) )
                % singulrity indicator passed but type not given.
                % Assume fractional poles or branches.
                if ( isSingEnd(1) )
                    % if singularity is at the left end point
                    obj.isSingEnd(1) = 1;
                    obj.singType{1} = 'branch';
                end
                
                if ( isSingEnd(2) )
                    % if singularity is at the right end point
                    obj.isSingEnd(2) = 1;
                    obj.singType{2} = 'branch';
                end
            end
            
            if ( nargin >= 4 && isempty(exponents) )
                % copy the information given about singularities in the current object
                obj.isSingEnd = isSingEnd;
                obj.singType = singType;                
            end          

            
            %%
            % Determine and factor out singular terms if exponents 
            % are not given
            if ( isempty(exponents) )
                obj.exponents = singfun.findSingExponents(op, obj.isSingEnd, obj.singType, pref);
                % update ISSINGEND and SINGTYPE based on EXPONENTS
                tol = pref.singfun.eps;
                if ( abs(obj.exponents(1)) < 100*tol )
                    % if the singularity exponent is below the tolerance level
                    % remove the singularity
                    obj.isSingEnd(1) = 0;
                    obj.singType{1} = 'none';
                end
                if ( abs(obj.exponents(2)) < 100*tol )
                    % if the singularity exponent is below the tolerance level
                    % remove the singularity
                    obj.isSingEnd(2) = 0;
                    obj.singType{2} = 'none';
                end
            end
            
            % update the operator based on the values in exponents.
            smoothOp = singfun.singOp2SmoothOp(op, obj.exponents, pref.singfun.eps);
            
            % Construct the smooth part of the SINGFUN object.
            % [TODO]: This will be replaced by the SMOOTHFUN constructor
            prefs = chebtech.pref('tech', 'cheb1', 'extrapolate', false);
            vscale = [];
            hscale = [];
            obj.smoothPart = chebtech.constructor(smoothOp, vscale, hscale, prefs);
        end
    end
    
    %% METHODS IMPLEMENTED BY THIS CLASS.
    methods
        
        % Complex conjugate of a SINGFUN.
        f = conj(f)
        
        % SINGFUN obects are not transposable.
        f = ctranspose(f)

        % Indefinite integral of a SINGFUN.
        f = cumsum(f, m, pref)

        % Derivative of a SINGFUN.
        f = diff(f, k, dim)
        
        % Extrapolate (for NaNs / Infs).
        [values, maskNaN, maskInf] = extrapolate(f)

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

        % Test if CHEBTECH objects are equal.
        out = isequal(f, g)

        % Test if a CHEBTECH is bounded.
        out = isfinite(f)

        % Test if a CHEBTECH is unbounded.
        out = isinf(f)

        % Test if a CHEBTECH has any NaN values.
        out = isnan(f)

        % True for real CHEBTECH.
        out = isreal(f)
        
        % True for zero CHEBTECH objects
        out = iszero(f)
        
        % Length of a CHEBTECH.
        len = length(f)

        % Convert a array-valued CHEBTECH into an ARRAY of CHEBTECH objects.
        g = mat2cell(f, M, N)

        % Global maximum of a CHEBTECH on [-1,1].
        [maxVal, maxPos] = max(f)

        % Global minimum of a CHEBTECH on [-1,1].
        [minVal, minPos] = min(f)

        % Global minimum and maximum on [-1,1].
        [vals, pos] = minandmax(f)

        % Subtraction of two CHEBTECH objects.
        f = minus(f, g)

        % Left matrix divide for CHEBTECH objects.
        X = mldivide(A, B)

        % Right matrix divide for a CHEBTECH.
        X = mrdivide(B, A)

        % Multiplication of CHEBTECH objects.
        f = mtimes(f, c)
        
        % Basic linear plot for CHEBTECH objects.
        varargout = plot(f, varargin)
        
        % Obtain data used for plotting a CHEBTECH object:
        data = plotData(f)

        % Addition of two CHEBTECH objects.
        f = plus(f, g)

        % Return the points used by a CHEBTECH.
        out = points(f)

        % QR factorisation of an array-valued CHEBTECH.
        [f, R, E] = qr(f, flag)

        % Right array divide for a CHEBTECH.
        f = rdivide(f, c, pref)

        % Real part of a SINGFUN.
        f = real(f)

        % Restrict a CHEBTECH to a subinterval.
        f = restrict(f, s)

        % Roots of a CHEBTECH in the interval [-1,1].
        out = roots(f, varargin)

        % Test an evaluation of the input OP against a CHEBTECH approx.
        pass = sampleTest(op, f)

        % Trim trailing Chebyshev coefficients of a CHEBTECH object.
        f = simplify(f, pref, force)

        % Size of a CHEBTECH.
        [siz1, siz2] = size(f, varargin)

        % Strict happiness check.
        [ishappy, epslevel, cutoff] = strictCheck(f, pref)

        % Definite integral of a CHEBTECH on the interval [-1,1].
        out = sum(f, dim)

        % CHEBTECH multiplication.
        f = times(f, g, varargin)
        
        % CHEBTECH obects are not transposable.
        f = transpose(f)

        % Unary minus of a CHEBTECH.
        f = uminus(f)

        % Unary plus of a CHEBTECH.
        f = uplus(f)

    end

    %% STATIC METHODS IMPLEMENTED BY THIS CLASS.
    methods ( Static = true )
        
        exponents = findSingExponents( op, isSingEnd, singType, pref )
        
        poleOrder = findPoleOrder( fvals, x, tol)
        
        barnchOrder = findBranchOrder( fvals, x, tol)
        
        op = singOp2SmoothOp( op, exponents, tol )
        
        % Retrieve and modify preferences for this class.
        prefs = pref(varargin)

    end

end
    
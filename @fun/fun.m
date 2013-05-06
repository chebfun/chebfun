classdef fun % (Abstract)
%FUN  Abstract FUN class for representing global functions on [a, b].

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FUN Class Description:
%
% The FUN class is an abstract class for representations of functions on the
% interval [a, b]. It acheives this my taking a onefun on [-1, 1] and applying
% a mapping.
%
% The current instances of FUNs are BNDFUNS and UNBNDFUNS. The former are used
% to represent functions on bounded domains, whereas the latter are able to
% represent some functions on unbounded domains.
%
% Class diagram: [chebfun] <>-- [<<FUN>>] <>-- [<<onefun>>]
%                                 ^   ^
%                                /     \
%                          [bndfun]   [unbndfun]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %% Properties of FUN objects.
    properties (Access = public)
        domain
        mapping
        onefun
    end
    
    %% CLASS CONSTRUCTOR:
    methods (Static)
        
        function obj = constructor(op, domain, hscale, vscale, pref)
            
            % Construct an empty fun:
            if ( nargin == 0 )
                op = [];
            end
            
            % Obtain preferences if none given:
            if ( nargin < 5 )
                pref = fun.pref;
            else
                pref = fun.pref(pref);
            end
           
            % Get domain if none given:
            if ( nargin < 2 )
                domain = pref.fun.domain;
            end

            % Get scales if none given:
            if ( nargin < 3 || isstruct(hscale) )
                if ( nargin > 2 && isstruct(hscale) )
                    pref = hscale; 
                end
                hscale = norm(domain, inf); 
                if ( isinf(hscale) )
                    hscale = 1; 
                end
            end
            if ( nargin < 4 )
                vscale = 0;
            end
            
            % Call constructor depending on domain:
            if ( ~any(isinf(domain)) )
                obj = bndfun(op, domain, hscale, vscale, pref);
            else
                obj = unbndfun(op, domain, hscale, vscale, pref);
            end
            
        end

    end
    
    %% STATIC METHODS IMPLEMENTED BY THIS CLASS.
    methods (Static = true)
        
        % Retrieve and modify preferences for this class.
        prefs = pref(varargin);
        
        % Edge detector.
        [edge, vscale] = detectedge(op, domain, hscale, vscale, pref, d)

    end
    
    %% ABSTRACT STATIC METHODS REQUIRED BY THIS CLASS.
    methods (Abstract = true, Static = true)

        % Map from [-1, 1] to the domain of the FUN.
        m = createMap(domain);  

    end
    
    %% METHODS IMPLEMENTED BY THIS CLASS.
    methods

%         % Convert an array of FUN objects into a array-valued FUN.
%         f = cell2mat(f)
        
%         % Plot (semilogy) the Chebyshev coefficients of a FUN object.
%         h = chebpolyplot(f, varargin)

        % Complex conjugate of a FUN.
        f = conj(f)
        
        % FUN objects are not transposable.
        f = ctranspose(f)

        % Evaluate a FUN.
        y = feval(f, x)

        % Flip columns of a vectorised FUN object.
        f = fliplr(f)
        
        % Imaginary part of a FUN.
        f = imag(f)

        % Compute the inner product of two FUN objects.
        out = innerProduct(f, g)

        % True for an empty FUN.
        out = isempty(f)

        % Test if FUN objects are equal.
        out = isequal(f, g)

        % Test if a FUN is bounded.
        out = isfinite(f)

        % Test if a FUN is unbounded.
        out = isinf(f)

        % Test if a FUN has any NaN values.
        out = isnan(f)

        % True for real FUN.
        out = isreal(f)
        
        % True for zero FUN objects
        out = iszero(f)
        
        % Length of a FUN.
        len = length(f)

        % Convert a array-valued FUN into an ARRAY of FUN objects.
        g = mat2cell(f, M, N)

        % Global maximum of a FUN on [a,b].
        [maxVal, maxPos] = max(f)

        % Global minimum of a FUN on [a,b].
        [minVal, minPos] = min(f)

        % Global minimum and maximum on [a,b].
        [vals, pos] = minandmax(f)

        % Subtraction of two FUN objects.
        f = minus(f, g)

        % Left matrix divide for FUN objects.
        X = mldivide(A, B)

        % Right matrix divide for a FUN.
        X = mrdivide(B, A)

        % Multiplication of FUN objects.
        f = mtimes(f, c)

        % Basic linear plot for FUN objects.
        varargout = plot(f, varargin)

        % Addition of two FUN objects.
        f = plus(f, g)

        % Polynomial coefficients of a FUN.
        out = poly(f)

        % QR factorisation of an array-valued FUN.
        [f, R, E] = qr(f, flag)

        % Right array divide for a FUN.
        f = rdivide(f, c, pref)

        % Real part of a FUN.
        f = real(f)

        % Roots of a FUN in the interval [a,b].
        out = roots(f, varargin)

%         % Trim trailing Chebyshev coefficients of a FUN object.
%         f = simplify(f, pref, force)

        % Size of a FUN.
        [siz1, siz2] = size(f, varargin)

        % FUN multiplication.
        f = times(f, g, varargin)
        
        % FUN obects are not transposable.
        f = transpose(f)

        % Unary minus of a FUN.
        f = uminus(f)

        % Unary plus of a FUN.
        f = uplus(f)

    end
end

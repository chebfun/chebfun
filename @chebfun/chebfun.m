classdef chebfun
%CHEBFUN    CHEBFUN class for representing piecewise smooth functions on [a, b].

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CHEBFUN Class Description:
%
% The CHEBFUN class is for representations of piecewise smooth functions on the
% interval [a, b].
%
% The CHEBFUN class is the main interface that a user will be using. Therefore,
% documentation should written appropriately, and checking of input arguments
% may be a little more forgiving.
%
% Class diagram: [jacfun] <>-- [CHEBFUN] <>-- [<<fun>>]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    properties (Access = public)
        domain
        funs
        impulses = [];
        isTransposed = 0;
    end

    methods

        function f = chebfun(varargin)
            % The main CHEBFUN constructor!

            % Return an empty CHEBFUN:
            if ( nargin == 0 || isempty(varargin{1}) )
                return
            end

            % Parse inputs:
            [op, dom, pref] = inputParser(varargin{:});

            if ( iscell(op) && all(cellfun(@(x) isa(x, 'fun'), op)) )
                % Construct a CHEBFUN from an array of FUN objects:

                if ( nargin > 1 )
                    error('CHEBFUN:chebfun:nargin', ...
                        'Only one input when passing an array of funs.')
                end
                f.funs = op;
                % Collect the domains together:
                dom = cellfun(@(f) get(f, 'domain'), f.funs, ...
                    'uniformOutput', false);
                f.domain = unique([dom{:}]);
                % Update values at jumps (first row of imps):
                f.impulses = chebfun.jumpVals(f.funs, f.domain, op);

            else
                % Construct from function_handle, numeric, or string input:

                % Convert string input to function_handle:
                if ( ischar(op) )
                    op = str2op(op);
                end

                % Can only return stats froma function_handle input:
                if ( ~isa(op, 'function_handle') )
                    pref.chebfun.stats = false;
                end
                % Initialise statistics and redefine operator:
                if ( pref.chebfun.stats )
                    op = @(x) statop(x, op);
                    op(NaN);
                end

                % Call the main constructor:
                [f.funs, f.domain, op] = chebfun.constructor(op, dom, pref);

                % Update values at jumps (first row of imps):
                f.impulses = chebfun.jumpVals(f.funs, f.domain, op);

                % Remove unnecessary breaks (but not those that were given):
                [ignored, index] = setdiff(f.domain, dom);
%                 f = merge(f, index', pref);

                % Reset statistics data:
                if ( pref.chebfun.stats )
                    op('foo');
                end

            end

        end

    end

    % Static methods implimented by CHEBFUN class.
    methods (Static = true)

        % Retrieve and modify preferences for this class.
        prefs = pref(varargin)

        % Splitting constructor.
        [funs, domain, op] = constructor(op, domain, pref)

        % Edge detector.
        [edge, vscale] = detectEdge(op, domain, hscale, vscale, derHandle)

        % Determine values of chebfun at breakpoints.
        vals = jumpVals(funs, ends, op)

        % Create a linear chebfun on a domain.
        out = x(dom, pref);

    end


    methods

        % Retrieve and modify preferences for this class.
        varargout = subsasgn(f, varargin);

    end

end

function op = str2op(op)
    % This is here as it's a clean function with no other variables hanging
    % around in the scope.
    depVar = symvar(op);
    if ( numel(depVar) ~= 1 )
        error('CHEBFUN:STR2OP:indepvars', ...
            'Incorrect number of independent variables in string input.');
    end
    op = eval(['@(' depVar{:} ')', op]);
end

function y = statop(x, op)
    % Allows statistics to be returns about construction from a function_handle.

    % Persistent variables:
    persistent store_x
    persistent caller

    % initialise them to empty cells:
    if ( isempty(store_x) )
        store_x = {};
        caller = {};
    end

    if ( isnan(x) )
        % Passing a NaN resets:
        store_x = {};
        caller = {};

    elseif ( ischar(x) )
        % Passing a string displays data:
        numEvals = sum(cellfun(@numel, store_x));
        fprintf('Total number of function evaluations = %d\n', numEvals);
        x = cellfun(@numel, store_x,'uniformoutput', false);
        disp([caller.' x.'])
        x = sort(vertcat(store_x{:}));
        hist(x, 9), shg
        store_x = {};
        caller = {};
    else
        if ( isempty(x) )
            return
        end
        % Else the operator is evaluated:
        y = op(x);
        %  and inputs are stored:
        store_x = [store_x x];
        tmp = dbstack;
        tmp(1:2) = [];
        if ( strncmpi(tmp(1).name, '@(x)', 4) )
            tmp(1) = [];
        end
        caller = [caller tmp(1).name];

    end

end

function [op, domain, pref] = inputParser(op, domain, varargin)
        % Parse inputs.

        args = varargin;
        if ( nargin == 1 )
            % chebfun(op)
            pref = chebfun.pref();
            domain = pref.chebfun.domain;
        elseif ( isstruct(domain) )
            % chebfun(op, pref)
            pref = domain;
            domain = pref.chebfun.domain;
        elseif ( nargin < 3 )
            % chebfun(op, domain)
            pref = chebfun.pref();
        elseif ( isstruct(varargin{1}) )
            % chebfun(op, domain, pref)
            pref = chebfun.pref(args{1});
            args(1) = [];
        elseif ( ~isnumeric(domain) )
            % chebfun(op, prop1, val1, ...)
            pref = chebfun.pref();
            args = [domain, args];
            domain = pref.chebfun.domain;
        else
            % chebfun(op, domain, prop1, val1, ...)
            pref = chebfun.pref;
        end

        % Obtain additional preferences:
        while ( ~isempty(args) )
            % Update these preferences:
            pref = chebfun.pref(pref, args{1}, args{2});
            args(1:2) = [];
        end

end



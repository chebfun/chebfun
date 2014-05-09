classdef cheboppref < chebpref
%CHEBOPPREF   Class for managing preferences for the Chebfun ODE suite.
%   CHEBOPPREF is a class for managing CHEBOP construction-time and solver
%   preferences, such as what solver is used for linear problem, the error or
%   residual tolerance for nonlinear problems, whether damped Newton iteration
%   should be performed for nonlinear problems, and how much information is to
%   be printed to the console during while the solver is active. 
%
% Available Preferences:
%
%   domain                     - Construction domain.
%     [-1, 1]
%
%     This sets the default domain that will be used for CHEBOP construction if
%     no domain argument is explicitly passed to the constructor.
%
%   discretization             - Discretization of linear problems
%     [@colloc2]
%     @ultraS
%
%     This options determines whether linear operators are discretized using
%     rectangular collocation methods or the ultraspherical method.
%
%   dimensionValues             - Increments in discretization sizes
%     [ 32    64   128   256   512   724  1024  1448]
%
%     This vector determines the number of gridpoints/coefficients used as
%     linear operators are discretized at finer and finer grids to resolve the
%     solution. For example, using the default value, a linear operator would
%     first be discretized at a 32 point grid, then a 64 point grid, up until a
%     1448 point grid.
%  
%   damped                      - Should Newton's method be damped?
%     [true]
%     false
%
%     If true, damped Newton iteration in function space is performed for
%     nonlinear problems. If false, undamped Newton iteration is performed, that
%     is, the solver will always take full Newton steps.
%
%   display                     - How much information is to be printed
%     'final'
%     ['iter']
%      'off'
%
%     If 'final', information is only printed after the solver of BVPs has
%     finished. If 'iter', information is printed at every Newton step. If
%     'off', no information is printed.
%
%   errTol                     - Error tolerance
%     [1e-10]
%
%     The termination criteria for the Newton iteration. The Newton iteration is
%     considered to have converged if the error estimate it computes is less
%     than the value of errTol.
%
%   lambdaMin                   - Minimum allowed step-size
%     [1e-6]
%
%     The value of lambdaMin determines the minimum allowed step-size that the
%     damped Newton iteration is allowed to take.
%
%   maxIter                     - Maximum number of Newton steps
%     [1e-6]
%
%   The maximum number of steps that the (damped) Newton iteration is allowed to
%   take, before it is considered to be non-convergent.
%
%   plotting                    - Plotting of intermediate Newton steps
%     DELAY
%     'on'
%     ['off']
%     'pause'
%
%   If plotting = 'on', the current iterate in the Newton solution is plotted at
%   every step, as well as the current Newton correction. If plotting = DELAY,
%   where DELAY has a numerical value, the iteration is paused and the plots are
%   shown for the time DELAY seconds. If plotting = 'pause', the iteration is
%   paused and the plots are shown until the user presses a button. If plotting
%   = 'off', no plots are shown during the Newton iteration.
%
% The default values for any of these preferences may be globally overridden
% using CHEBOPPREF.SETDEFAULTS(); see the documentation for that function for
% further details.
%
% Constructor inputs:
%   P = CHEBOPPREF() creates a CHEBOPPREF object with the default values of the
%   preferences.  For a list of all available preferences, see above.
%
%   P = CHEBOPPREF(Q), where Q is a MATLAB structure uses the field/value pairs
%   in Q to set the properties of P.  If a field of Q has a name which matches
%   a property of P, the value of that property of P is set to the value
%   associated to that field in Q.  If a field of Q has a name that does not
%   correspond to a known preference, then an error is thrown.
%
%   P = CHEBOPPREF(Q), where Q is a CHEBOPPREF, sets P to be a copy of Q.
%
% See also CHEBFUNPREF.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% TODO:  Further documentation of CHEBOPPREF preferences.

    methods

        function outPref = cheboppref(inPref)
            if ( (nargin == 1) && isa(inPref, 'cheboppref') )
                outPref = inPref;
                return
            elseif ( nargin < 1 )
                inPref = struct();
            end

            % Initialize default preference values.
            outPref.prefList = cheboppref.manageDefaultPrefs('get');

            % Copy fields from q, merging incomplete substructures.
            for field = fieldnames(inPref).'
                if ( isfield(outPref.prefList, field{1}) )
                    if ( isstruct(outPref.prefList.(field{1})) )
                        outPref.prefList.(field{1}) = ...
                            chebpref.mergePrefs(outPref.prefList.(field{1}), ...
                            inPref.(field{1}));
                    else
                        outPref.prefList.(field{1}) = inPref.(field{1});
                    end
                else
                    error('CHEBOPPREF:cheboppref:badPref', ...
                        'Unrecognized preference name.');
                end
            end
        end

       function display(pref)
       %DISPLAY   Display a CHEBOPPREF object.
       %   DISPLAY(PREF) prints out a list of the preferences stored in the
       %   CHEBOPPREF object PREF.

            % Compute the screen column in which pref values start.
            valueCol = 24; % length('    enableSingularityDetection:   ');

            % A subfunction to pad strings for formatting.
            function s = padString(s)
            %PADSTRING   Add whitespace to string for formatting.
                s = [s repmat(' ', 1, valueCol - length(s))];
            end

            % Print values of "known" preferences.
            prefList = pref.prefList;

            fprintf('cheboppref object with the following preferences:\n');
            fprintf([padString('    domain:') '[%g, %g]\n'], ...
                prefList.domain(1), prefList.domain(end));
            fprintf([padString('    discretization:') '%s\n'], ...
                func2str(prefList.discretization));
            fprintf([padString('    dimensionValues:') '%s\n'], ...
                num2str(prefList.dimensionValues));
            fprintf([padString('    damped:') '%d\n'], ...
                prefList.damped);
            fprintf([padString('    display:') '%s\n'], ...
                prefList.display);
            fprintf([padString('    errTol:') '%g\n'], ...
                prefList.errTol);
            fprintf([padString('    lambdaMin:') '%g\n'], ...
                prefList.lambdaMin);
            fprintf([padString('    maxIter:') '%d\n'], ...
                prefList.maxIter);
            fprintf([padString('    plotting:') '%s\n'], ...
                prefList.plotting);
        end

    end

    methods ( Static = true )
        function pref = getFactoryDefaults(getFactory)
        %GETFACTORYDEFAULTS   Get factory default preferences.
        %   PREF = CHEBOPPREF.GETFACTORYDEFAULTS() returns a CHEBOPPREF
        %   object with the preferences set to their factory defaults,
        %   irrespective of the currently defined values of the default
        %   preferences.  This function is useful if the user wishes to
        %   solve ODEs with CHEBOP using the factory defaults when other
        %   user-set defaults are currently in force.
        %
        % See also GETDEFAULTS, SETDEFAULTS.

            fd = cheboppref.factoryDefaultPrefs();
            pref = cheboppref(fd);
        end

        function pref = getDefaults()
        %GETDEFAULTS   Get default preferences.
        %   PREF = CHEBOPPREF.GETDEFAULTS() returns a CHEBOPPREF object with
        %   the preferences set to the currently stored default values.  It is
        %   equivalent to PREF = CHEBOPPREF().
        %
        % See also GETFACTORYDEFAULTS, SETDEFAULTS.

            pref = cheboppref();
        end

        function setDefaults(varargin)
        %SETDEFAULTS   Set default preferences.
        %   CHEBOPPREF.SETDEFAULTS(PREF1, VAL1, PREF2, VAL2, ...) sets the
        %   default values for the preferences whose names are stored in the
        %   strings PREF1, PREF2, ..., etc. to VAL1, VAL2, ..., etc.  All
        %   subsequently constructed CHEBOPPREF objects will use these values
        %   as the defaults.
        %
        %   CHEBOPPREF.SETDEFAULTS(PREF) sets the default values to the
        %   preferences stored in the CHEBOPPREF object PREF.  PREF can also
        %   be a MATLAB structure, in which case it is converted to a
        %   CHEBOPPREF as described in the documentation for the CHEBOPPREF
        %   constructor first.
        %
        %   CHEBOPPREF.SETDEFAULTS('factory') resets the default preferences to
        %   their factory values.
        %
        % See also GETDEFAULTS, GETFACTORYDEFAULTS.

            % The reason we don't just use manageDefaults as the second
            % argument to chebpref.setDefaults and wrap it in an additional
            % anonymous function instead is to get around what seems to be a
            % bug in MATLAB.  See commit messages for more information.
            manageDefaults = @cheboppref.manageDefaultPrefs;
            chebpref.setDefaults(@(inPref) cheboppref(inPref), ...
                @(varargin) manageDefaults(varargin{:}), varargin{:});
        end
    end

    methods ( Static = true, Access = private )

        function varargout = manageDefaultPrefs(varargin)
        %MANAGEDEFAULTPREFS   Private method for handling default preferences.
        %   CHEBOPPREF.MANAGEDEFAULTPREFS('get') returns a structure suitable
        %   for storing in the prefList property of a CHEBOPPREF with all of
        %   the currently stored default preferences suitable for initializing
        %   a CHEBOPPREF object.
        %
        %   CHEBOPPREF.MANAGEDEFAULTPREFS('set-factory') restores the default
        %   preferences to their "factory" values.
        %
        %   CHEBOPPREF.MANAGEDEFAULTPREFS('set', PREFLIST) sets the default
        %   values to those stored in the structure PREFLIST.  PREFLIST should
        %   be a structure suitable for use as a CHEBOPPREF prefList.
        %
        %   CHEBOPPREF.MANAGEDEFAULTPREFS('set', PREF1, VAL1, PREF2, VAL2, ...)
        %   sets the default values for PREF1, PREF2, ..., etc. to VAL1, VAL2,
        %   ..., etc.

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Developer notes:
        %  - MATLAB has no equivalent to what might be called a "static" class
        %    variable in other languages, so a persistent variable is the best
        %    we can do for providing this feature.  Persistent variables are
        %    local to a specific function, so we can have only a single
        %    function for managing it.  As a result, this function has a mildly
        %    awkward syntax and so is not user-facing.
        %  - More importantly, this function is also not user-facing because
        %    its inputs and outputs depend on the internal representation of a
        %    CHEBOPPREF as a MATLAB structure, and that's not something with
        %    which anyone outside of CHEBOPPREF should be concerned.
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            persistent defaultPrefs;

            if ( isempty(defaultPrefs) )
                defaultPrefs = cheboppref.factoryDefaultPrefs();
            end

            if ( strcmp(varargin{1}, 'get') )
                varargout{1} = defaultPrefs;
            elseif ( strcmp(varargin{1}, 'set-factory') )
                defaultPrefs = cheboppref.factoryDefaultPrefs();
            elseif ( strcmp(varargin{1}, 'set') )
                    varargin(1) = [];
                if ( isstruct(varargin{1}) )
                    defaultPrefs = varargin{1};
                else
                    while ( ~isempty(varargin) )
                        prefName = varargin{1};
                        prefValue = varargin{2};
                        if ( isfield(defaultPrefs, prefName) )
                            defaultPrefs.(prefName) = prefValue;
                        else
                            error('CHEBOPPREF:manageDefaultPrefs:badPref', ...
                                'Unrecognized preference name.');
                        end
                        varargin(1:2) = [];
                    end
                end
            end
        end

        function factoryPrefs = factoryDefaultPrefs()
        %FACTORYDEFAULTPREFS   Get structure of factory default preferences.
        %   S = CHEBOPPREF.FACTORYDEFAULTPREFS() returns a structure suitable
        %   for storing in the prefList property of a CHEBOPPREF object that
        %   contains all of the "factory default" values of the CHEBOP
        %   preferences.

            factoryPrefs.domain = [-1 1];
            factoryPrefs.discretization = @colloc2;
            factoryPrefs.scale = NaN;
            factoryPrefs.dimensionValues = [32 64 128 256 512 724 1024 1448 2048];
            factoryPrefs.damped = 1;
            factoryPrefs.display = 'off';
            factoryPrefs.errTol = 1e-10;
            factoryPrefs.lambdaMin = 1e-6;
            factoryPrefs.maxIter = 25;
            factoryPrefs.plotting = 'off';
        end

    end

end

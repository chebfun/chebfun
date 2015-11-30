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
%     ['values']
%     'coeffs'
%     @chebcolloc1
%     @chebcolloc2
%     @ultraS
%     @trigcolloc
%     @trigspec
%
%     This options determines whether linear operators are discretized using
%     rectangular collocation methods or the ultraspherical method. Please
%     observe that
%         * 'values' and 'coeffs' are convenient ways of specifying the 
%           @chebcolloc2 and @ultraS options respectively (when the boundary 
%           conditions are not periodic), and @trigcolloc and @trigspec 
%           respectively (when the boundary conditions are periodic).
%         * The @trigcolloc/@trigspec options are only supported for problems
%           that are specified to have periodic boundary conditions. 
%         * Specifying the @chebcolloc1 option causes the CHEBFUN solution
%           returned to be based on the @chebtech1 tech. The @chebtech2/@ultraS
%           option causes the CHEBFUN solution returned to be based on the
%           @chebtech2 tech. The @trigcolloc/@trigspec option causes the 
%           CHEBFUN solution to be periodic, based on the @trigtech tech.
%        
%   damping                     - Should Newton's method be damped?
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
%   errTol                      - Error tolerance
%     [1e-10]
%
%     The termination criteria for the Newton iteration. The Newton iteration is
%     considered to have converged if the error estimate it computes is less
%     than the value of errTol.
%
%   happinessCheck              - Routine for checking that solution converged
%     [@standardCheck]
%     @basicCheck
%     @plateauCheck
%     @classicCheck
%     @looseCheck
%     @strictCheck
%     @happinessCheck
%     @linopV4Check
%
%     This options determines which routine is used to determine that the
%     approximate solution has converged. Any of the above options may be
%     used, as well as any user defined function handle that conforms to 
%     the happinessCheck standards.
%
%   ivpAbsTol                    - Absolute tolerance for the ivpSolver
%     [1e5*eps]
%
%     This options specifies the option for the absolute tolerance passed as an
%     option to the built-in MATLAB ODE solver when solving IVPs.
%
%   ivpRelTol                    - Relavtive tolerance for the ivpSolver
%     [100*eps]
%
%     This options specifies the option for the relative tolerance passed as an
%     option to the built-in MATLAB ODE solver when solving IVPs.
%
%   ivpRestartSolver             - Restart IVP solvers at breakpoints
%     false
%     [true]
%
%     This option specifies whether the MATLAB built in solvers should be
%     restarted at breakpoints. That is, whether each subinterval of a piecewise
%     problem will get integrated separately. This can be very useful for e.g.
%     short forcing pulses, which otherwise might get overlooked.
%
%   ivpSolver                  - Solver for IVPs
%     ['ode113']
%     'ode15s'
%     'ode45'
%     'values'
%     'coeffs'
%
%     This options determines which of the MATLAB built-in IVP solvers is used
%     for solving IVPs posed with the CHEBOP class. Any option of
%     CHEBOPPREF.discretization (see above) is allowed, which causes IVPs to be
%     solved globally via spectral methods, rather than reformulating them as
%     first-order problems and then solved via time-stepping method.
%
%   lambdaMin                   - Minimum allowed step-size for Newton's method
%     [1e-6]
%
%     The value of lambdaMin determines the minimum allowed step-size that the
%     damped Newton iteration is allowed to take.
%
%   maxDimension
%     [4096]
%
%     The maximum number of gridpoints/coefficients used as linear operators are
%     discretized at finer and finer grids to resolve the solution. The
%     intermediate values for the discretization between cheboppref.minDimension
%     and cheboppref.maxDimension depend on the discretization used for the
%     operator.
%
%   maxIter                     - Maximum number of Newton steps
%     25
%
%   The maximum number of steps that the (damped) Newton iteration is allowed to
%   take, before it is considered to be non-convergent.
%
%   minDimension
%     [32]
%
%     The minimum number of gridpoints/coefficients used as linear operators are
%     discretized at finer and finer grids to resolve the solution. The
%     intermediate values for the discretization between cheboppref.minDimension
%     and cheboppref.maxDimension depend on the discretization used for the
%     operator.
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
%   vectorize                   - Automatic vectorization of anon. functions
%     [true]
%     false
%
%   Determines whether the CHEBOP class should try to automatically try to
%   vectorize anonymous functions used for describing the differential equation
%   and boundary condition(s).
%
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

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% TODO:  Further documentation of CHEBOPPREF preferences.

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CLASS CONSTRUCTOR:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods ( Access = public, Static = false )

        function outPref = cheboppref(inPref, varargin)
            if ( (nargin == 1) && isa(inPref, 'cheboppref') )
                outPref = inPref;
                return
            elseif ( nargin < 1 )
                inPref = struct();
            elseif ( ischar(inPref) )
                if ( nargin == 1 )
                    error('CHEBFUN:CHEBOPPREF:cheboppref:deprecated', ...
                        ['cheboppref() no longer supports queries of ', ...
                         'the form cheboppref(''prop'').\n', ...
                         'Please use cheboppref().prop.']);
                else
                    error('CHEBFUN:CHEBOPPREF:cheboppref:deprecated', ...
                        ['cheboppref() no longer supports assignment ', ...
                         'via cheboppref(''prop'', val).\n', ...
                         'Please use cheboppref.setDefaults(''prop'', val).']);
                end
            elseif ( nargin > 1 )
                error('CHEBFUN:CHEBOPPREF:cheboppref:inputs', ...
                    'Too many input arguments.')
            end

            % Initialize default preference values.
            outPref.prefList = cheboppref.manageDefaultPrefs('get');

            % Copy fields from q, merging incomplete substructures.
            for field = fieldnames(inPref).'
                field1 = field{1};
                if ( isfield(outPref.prefList, field1) )
                    if ( isstruct(outPref.prefList.(field1)) )
                        outPref.prefList.(field1) = ...
                            chebpref.mergePrefStructs(...
                                outPref.prefList.(field1), ...
                                inPref.(field1));
                    else
                        outPref.prefList.(field1) = inPref.(field1);
                    end
                else
                    error('CHEBFUN:CHEBOPPREF:cheboppref:badPref', ...
                        'Unrecognized preference name.');
                end
            end
        end
        
    end
    

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CLASS METHODS:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods ( Access = public, Static = false )

       function display(pref)
       %DISPLAY   Display a CHEBOPPREF object.
       %   DISPLAY(PREF) prints out a list of the preferences stored in the
       %   CHEBOPPREF object PREF.

            % Compute the screen column in which pref values start.
            valueCol = 24; % length('    blowup:   ');

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
            if ( isa(prefList.discretization,'function_handle') )
                fprintf([padString('    discretization:') '%s\n'], ...
                    func2str(prefList.discretization));
            elseif ( isa(prefList.discretization,'char') )
                fprintf([padString('    discretization:') '%s\n'], ...
                    prefList.discretization);
            end
  
            fprintf([padString('    damping:') '%d\n'], ...
                prefList.damping);
            fprintf([padString('    display:') '%s\n'], ...
                prefList.display);
            fprintf([padString('    errTol:') '%g\n'], ...
                prefList.errTol);
            fprintf([padString('    happinessCheck:') '%s\n'], ...
                func2str(prefList.happinessCheck));
            fprintf([padString('    ivpAbsTol:') '%g\n'], ...
                prefList.ivpAbsTol);
            fprintf([padString('    ivpRelTol:') '%g\n'], ...
                prefList.ivpRelTol);
            fprintf([padString('    ivpRestartSolver:') '%d\n'], ...
                prefList.ivpRestartSolver);
            fprintf([padString('    ivpSolver:') '%s\n'], ...
                func2str(prefList.ivpSolver));
            fprintf([padString('    lambdaMin:') '%g\n'], ...
                prefList.lambdaMin);
            fprintf([padString('    maxDimension:') '%d\n'], ...
                prefList.maxDimension);
            fprintf([padString('    maxIter:') '%d\n'], ...
                prefList.maxIter);
            fprintf([padString('    minDimension:') '%d\n'], ...
                prefList.minDimension);
            fprintf([padString('    plotting:') '%s\n'], ...
                prefList.plotting);
            fprintf([padString('    vectorize:') '%i\n'], ...
                prefList.vectorize);
       end

        function pref = subsasgn(pref, ind, val)
        %SUBSASGN   Subscripted assignment for CHEBOPPREF.
        %   P.PROP = VAL, where P is a CHEBOPPREF object, assigns the value
        %   VAL to the CHEBOPPREF property PROP stored in P.  If PROP is not a
        %   CHEBOPPREF property, an error will be thrown.
        %
        %   CHEBOPPREF does not support any other subscripted assignment types,
        %   including '()' and '{}'.
            
            % Support user-friendlier syntax for specifying discretization
            % choice:
            val = cheboppref.parseDiscretization(val);
            
            % Support user-friendlier syntax for specifying IVP solver choice:
            val = cheboppref.parseIVPsolver(val);
            
            % Call the superclass method.
            pref = subsasgn@chebpref(pref, ind, val);
        end 
       
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% STATIC METHODS:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    methods ( Access = public, Static = true )
        
        function pref = getFactoryDefaults()
        %GETFACTORYDEFAULTS   Get factory default preferences.
        %   PREF = CHEBOPPREF.GETFACTORYDEFAULTS() returns a CHEBOPPREF
        %   object with the preferences set to their factory defaults,
        %   irrespective of the currently defined values of the default
        %   preferences.  This function is useful if the user wishes to
        %   solve ODEs with CHEBOP using the factory defaults when other
        %   user-set defaults are currently in force.
        %
        % See also SETDEFAULTS.

            fd = cheboppref.factoryDefaultPrefs();
            pref = cheboppref(fd);
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
        % See also GETFACTORYDEFAULTS.

            % The reason we don't just use manageDefaults as the second
            % argument to chebpref.setDefaults and wrap it in an additional
            % anonymous function instead is to get around what seems to be a
            % bug in MATLAB.  See commit messages for more information.
            manageDefaults = @cheboppref.manageDefaultPrefs;
            chebpref.setDefaults(@(inPref) cheboppref(inPref), ...
                @(varargin) manageDefaults(varargin{:}), varargin{:});
        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% PRIVATE STATIC METHODS
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
        % DEVELOPER NOTE:
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
                        
                        % Support user-friendlier syntax for specifying
                        % discretization choice:
                        prefValue = cheboppref.parseDiscretization(prefValue);
                        prefValue = cheboppref.parseHappinessCheck(prefValue);
                        prefValue = cheboppref.parseIVPsolver(prefValue);
                        if ( isfield(defaultPrefs, prefName) )
                            defaultPrefs.(prefName) = prefValue;
                        else
                            error('CHEBFUN:CHEBOPPREF:cheboppref:badPref', ...
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
            factoryPrefs.discretization = 'values';
            factoryPrefs.scale = NaN;
            factoryPrefs.damping = 1;
            factoryPrefs.display = 'off';
            factoryPrefs.errTol = 1e-10;
            factoryPrefs.happinessCheck = @standardCheck;
            factoryPrefs.ivpAbsTol = 1e5*eps;
            factoryPrefs.ivpRelTol = 100*eps;
            factoryPrefs.ivpRestartSolver = true;
            factoryPrefs.ivpSolver = @chebfun.ode113;
            factoryPrefs.lambdaMin = 1e-6;
            factoryPrefs.maxDimension = 4096;
            factoryPrefs.maxIter = 25;
            factoryPrefs.minDimension = 32;
            factoryPrefs.plotting = 'off';
            factoryPrefs.vectorize = true;
        end
        
        function val = parseDiscretization(val)
        %PARSEDISCRETIZATION    Allow different syntax for specifying
        %                       discretization.
            
            % We want to allow user-friendly syntax for specifying the
            % discretization (#433). So check whether we have some of the
            % strings we want to allow, and convert them to the correct function
            % handle:
            if ( any(strcmpi(val, {'ultraspherical', 'ultraS'})) )
                warning('CHEBOPPREF:PARSEDISCRETIZATION', ...
                    ['''ULTRAS''/''ULTRASPHERICAL'' is deprecated. \n' ...
                    'Please use ''COEFFS''/@ultraS.']);
                val = @ultraS;
                
            elseif ( any(strcmpi(val, {'chebcolloc2', 'collocation', 'colloc2'})) )
                warning('CHEBOPPREF:PARSEDISCRETIZATION', ...
                    ['''COLLOCATION''/''COLLOC2''/''CHEBCOLLOC2'' is deprecated. \n' ...
                    'Please use ''VALUES''/@chebcolloc2.']);
                val = @chebcolloc2;
                
            elseif ( any(strcmpi(val, {'chebcolloc1', 'colloc1'})) )
                warning('CHEBOPPREF:PARSEDISCRETIZATION', ...
                    ['''COLLOC1''/''CHEBCOLLOC1'' is deprecated. \n' ...
                    'Please use ''VALUES''/@chebcolloc2.']);
                val = @chebcolloc1;
                
            elseif ( any(strcmpi(val, {'trigcolloc', 'periodic'})) )
                warning('CHEBOPPREF:PARSEDISCRETIZATION', ...
                    ['''TRIGCOLLOC''/''PERIODIC'' is deprecated. \n' ...
                    'Please use ''VALUES''/@trigcolloc.']);
                val = @trigcolloc;
            end
                
        end
        
        function val = parseHappinessCheck(val)
        %PARSEHAPPINESSCHECK    Allow different syntax for specifying
        %                       happinessCheck.
            
            % handle:
            if ( any(strcmpi(val, {'classic', 'classicCheck'})) )
                val = @classicCheck;
                
            elseif ( any(strcmpi(val, {'plateau', 'plateauCheck'})) )
                val = @plateauCheck;
                
            elseif ( any(strcmpi(val, {'strict', 'strictCheck'})) )
                val = @strictCheck;
                 
            elseif ( any(strcmpi(val, {'loose', 'looseCheck'})) )
                val = @looseCheck;
                 
            elseif ( any(strcmpi(val, {'happiness', 'happinessCheck'})) )
                val = @happinessCheck;
                 
            elseif ( any(strcmpi(val, {'linopV4', 'linopV4Check'})) )
                val = @linopV4Check;
                 
            end
                
        end
        
        function val = parseIVPsolver(val)
        %PARSEIVPSOLVER   Allow different syntax for specifying the IVPsolver.
            
            % Check whether we got pref.ivpSolver = @ode113/@ode45/@ode15s, that
            % is, a function handle, but not the CHEBFUN overload of it.
            if ( isa(val, 'function_handle') && ...
                    any(strcmpi(func2str(val), {'ode113', 'ode15s', 'ode45'})) )
                val = eval(['@chebfun.', func2str(val)]);
                
            % Check whether we got a string argument, e.g. 
            % pref.ivpSolver = 'ode113'.
            elseif ( any(strcmpi(val, {'ode113', 'ode15s', 'ode45'})) )
                val = eval(['@chebfun.', val]);
            end
        end

    end
    
end

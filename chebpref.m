classdef chebpref
%CHEBPREF   Class for managing Chebfun preferences.
%   CHEBPREF is a class for managing Chebfun construction-time preferences such
%   as the construction tolerance, whether or not to perform breakpoint and
%   singularity detection, and the various options that those features require.
%   These objects can be supplied to the CHEBFUN constructor (as well as the
%   constructors of other classes in the CHEBFUN system), which will interpret
%   them and adjust the construction process accordingly.
%
% Available Preferences:
%
%   domain                     - Construction domain.
%    [-1, 1]
%
%      This sets the default domain that will be used for CHEBFUN and/or FUN
%      construction if no domain argument is explicitly passed to the
%      constructor.
%
%   enableBreakpointDetection  - Enable/disable breakpoint detection.
%     true
%    [false]
%
%     If true, breakpoints between FUNS may be introduced where a discontinuity
%     in a function or a low-order derivative is detected or if a global
%     representation will be too long.  If false, breakpoints will be
%     introduced only at points where discontinuities are being created (e.g.,
%     by ABS(F) at points where a CHEBFUN F passes through zero).
%
%   breakpointPrefs            - Preferences for breakpoint detection.
%
%      splitMaxLength          - Maximum FUN length.
%       [129]
%
%         This is the maximum length of a single FUN (e.g., the number of
%         Chebyshev points used for FUNs based on Chebyshev polynomial
%         interpolation) allowed by the constructor when breakpoint detection
%         is enabled.
%
%      splitMaxTotalLength     - Maximum total CHEBFUN length.
%       [6000]
%
%         This is the maximum total length of the CHEBFUN (i.e., the sum of the
%         lengths of all the FUNs) allowed by the constructor when breakpoint
%         detection is enabled.
%
%   enableSingularityDetection - Enable/disable singularity detection.
%     true
%    [false]
%
%      If true, the constructor will attempt to detect and factor out
%      singularities, (e.g., points where a function or its derivatives become
%      unbounded). If false, breakpoints will be introduced only at points where
%      singularities are being created, (e.g., by SQRT(F) at points where a
%      CHEBFUN F passes through zero). See SINGFUN for more information.
%
%   enableDeltaFunctions - Enable delta functions.
%     true
%    [false]
%       
%      If true, the deltafun class will be invoked to manage any delta 
%      functions present in the object.
%   deltaPrefs                 - Preferences for delta functions.
%
%      deltaTol                - Tolerance for magnitude of delta functions.
%       [1e-11]
%
%         This is the tolerance up to which delta functions will be negligible.
%
%      proxomityTol          - Minimum distance between delta functions. 
%       [1e-11]
%         If two delta functions are located closer than this tolerance, they 
%         will be merged.
%
%   scale                      - The vertical scale constructor should use.
%    [0]
%
%      Typically the CHEBFUN constructor will resolve relative to a vertical
%      scale determined by it's own function evaluations. However, in some
%      situations one would like to the resolve relative to a fixed vertical
%      scale. This can be set using this preference.
%
%   singPrefs                  - Preferences for singularity detection.
%
%      exponentTol             - Tolerance for exponents.
%       [1.1*1e-11]
%
%         This is the tolerance up to which the detector will try to resolve
%         the singularity exponents.
%
%      maxPoleOrder            - Maximum pole order.
%       [20]
%
%         Maximum order of the pole that the singularity detector can find.
%
%      exponents               - Exponents at the end points.
%       [[ ]]
%
%         If exponents are supplied by the user from CHEBFUN or FUN levels, the
%         default value (empty) is replaced and passed to singfun constructor.
%         If no information about the singularity is specified by the user, then
%         the singularity detection is triggered to find the exact pole order
%         (which is integer) or singularity order (which is fractional) when
%         enableSingularityDetection is set TRUE. When exponents are given by
%         the user, there is no need to specify the 'singType' field, as any
%         information in singType will be ignored by the SINGFUN constructor.
%         For a piecewise smooth CHEBFUN, the number of exponents should be
%         given in pair with each pair corresponds to the ends of a piece. For
%         example,
%
%         dom = [-2 -1 0 1];
%         op1 = @(x) sin(x);
%         op2 = @(x) 1./(1+x);
%         op3 = @(x) x+1;
%         op = {op1, op2, op3};
%         pref = chebpref();
%         pref.singPrefs.exponents = [0 0 -1 0 0 0];
%         f = chebfun(op, dom, pref);
%
%         Note that syntax in Chebfun v4 is still supported. So the example
%         above can be exercised as below:
%
%         dom = [-2 -1 0 1];
%         op1 = @(x) sin(x);
%         op2 = @(x) 1./(1+x);
%         op3 = @(x) x+1;
%         op = {op1, op2, op3};
%         f = chebfun(op, dom, 'exps', [0 0 -1 0 0 0]);
%
%         For the cases where the CHEBFUN has more than one piece, if the size
%         of the given exponents is 1x2, then the CHEBFUN constructor will take
%         them as the exponent for the left endpoint of the first piece and the
%         and that for the right endpoint of the last piece. The exponents for
%         all other interior endpoints are simply assumed zeros. For example,
%
%         pref = chebpref();
%         pref.singPrefs.exponents = [-1 0];
%         pref.enableBreakpointDetection = 1;
%         f = chebfun(@(x) sin(100*x)./(x+2), [-2 7], pref)
%
%         The equivalent syntax in Chebfun v4 fashion is still valid:
%
%         f = chebfun(@(x) sin(100*x)./(x+2), [-2 7], 'splitting', 'on', ...
%             'exps', [-1 0])
%
%      singType                - Type of singularities.
%       [{ }]
%
%         The information provided in singType helps the singularity detector to
%         determine the order of the singularities more efficiently and save
%         some construction time. If the default, i.e. an empty cell, is
%         replaced by a user-supplied 2*N cell with entries being any of 'none',
%         'pole', 'sing', and 'root' where N is the number of smooth pieces,
%         i.e. FUNS, then this cell is passed to the SINGFUN constructor to
%         speed up the singularity detection. Here, 'none', 'pole', 'sing', and
%         'root' correspond to free of any kind of singularities, integer pole,
%         fractional singularity, and root at the end point with order less than
%         1, respectively. With the default empty cell, the SINGFUN constructor
%         will assume fractional singularities. For instance, setting singType
%         to {'pole', 'sing'} tells the singularity detector to search for poles
%         at the left endpoint of an interval and arbitrary singularities at the
%         right endpoint. For example,
%
%         dom = [-1 1];
%         op = @(x) (x - dom(1)).^-0.5.*sin(x);
%         pref = chebpref();
%         pref.singPrefs.singType = {'sing', 'none'};
%         f = chebfun(op, dom, pref);
%
%         Syntactically, chebfun constructor supports automatic singularity
%         detection for piecewise smooth CHEBFUN. That is, the users can specify
%         a series of the strings described above in pairs with each pair
%         corresponds to the endpoints of a subinterval. For example, if one
%         want to construct a CHEBFUN definied in [-1 0 1] with poles at -1 and
%         1 and fractional singularity on each side of 0, then the series of
%         string passed to the CHEBFUN constructor should be {'pole', 'sing',
%         'sing', 'pole'}. With these information, the CHEBFUN constructor and
%         consequently the SINGFUN will try to find the exact order of the
%         singularities. However, the SINGFUN constructor may not succeed for
%         most of the cases due to the unsatisfactory performance of the
%         current singularity detector.
%
%   scale                      - The vertical scale the constructor should use.
%    [0]
%
%      Typically the CHEBFUN constructor will resolve relative to a vertical
%      scale determined by it's own function evaluations. However, in some
%      situations one would like to the resolve relative to a fixed vertical
%      scale. This can be set using this preference.
%
%   tech                       - Representation technology.
%    ['chebtech']
%
%      Sets the underlying representation technology used to construct the FUNs.
%
%   techPrefs                  - Preferences for the tech constructor.
%
%      This is a structure of preferences that will be passed to the constructor
%      for the underlying representation technology.  See, for example,
%      CHEBTECH.TECHPREF for preferences accepted by the default CHEBTECH
%      technology.  Additionally, all techs are required to accept the following
%      preferences:
%
%      eps                     - Construction tolerance.
%       [2^(-52)]
%
%        Specifies the relative tolerance to which the representation should be
%        constructed.
%
%      maxLength               - Maximum representation length.
%       [65537]
%
%        Maximum length of the underlying representation.
%
%      exactLength             - Exact representation length.
%       [NaN]
%
%        Exact length of the underlying representation to be used.  A NaN value
%        indicates that the tech is free to choose the length (up to maxLength),
%        e.g., as the basis of an adaptive construction procedure.
%
%      extrapolate             - Extrapolate endpoint values.
%        true
%       [false]
%
%        If true, the tech should avoid direct evaluation of the function at
%        the interval endpoints and "extrapolate" the values at those points if
%        needed.  It should also extrapolate the values of any points for which
%        the function being sampled returns NaN.
%
%      sampleTest              - Test accuracy at arbitrary point.
%       [true]
%        false
%
%        If true, the tech should check an arbitrary point for accuracy to
%        ensure that behavior hasn't been missed, e.g., due to undersampling.
%
% The default values for any of these preferences may be globally overridden
% using CHEBPREF.SETDEFAULTS(); see the documentation for that function for
% further details.
%
% Constructor inputs:
%   P = CHEBPREF() creates a CHEBPREF object with the default values of the
%   preferences.  For a list of all available preferences, see above.
%
%   P = CHEBPREF(Q), where Q is a MATLAB structure uses the field/value pairs
%   in Q to set the properties of P.  If a field of Q has a name which matches
%   a property of P, the value of that property of P is set to the value
%   associated to that field in Q.  Any fields of Q that are not properties of
%   P are interpreted as preferences for the constructor of the underlying
%   representation technology and are placed in P.TECHPREFS.  The exceptions to
%   this are the fields BREAKPOINTPREFS, SINGPREFS, and TECHPREFS.  If Q has
%   fields with these names, they will be assumed to be MATLAB structures and
%   will be "merged" with the structures of default preferences stored in the
%   properties of the same names in P using CHEBPREF.MERGEPREFS().
%
%   P = CHEBPREF(Q), where Q is a CHEBPREF, sets P to be a copy of Q.
%
% Notes:
%   When building a CHEBPREF from a structure using the second calling syntax
%   above, one should take care to ensure that preferences for the underlying
%   representation technology are specified once and only once; e.g., do not
%   simultaneously set Q.MYPREF = 1 and Q.TECHPREFS.MYPREF = 2.  The value of
%   P.TECHPREFS.MYPREF that gets set from P = CHEBPREF(Q) in this circumstance
%   is implementation-defined.
%
% Examples:
%   Create a CHEBPREF for building a CHEBFUN based on CHEBTECH (default) with
%   breakpoint detection, a splitting length of 257 (pieces of polynomial degree
%   256, and a custom CHEBTECH refinement function:
%      p.enableBreakpointDetection = true;
%      p.breakpointPrefs.splitLength = 257;
%      p.techPrefs.refinementFunction = @custom;
%      pref = chebpref(p);
%
%   Same thing with a slightly shorter syntax:
%      p.enableBreakpointDetection = true;
%      p.breakpointPrefs.splitLength = 257;
%      p.refinementFunction = @custom;
%      pref = chebpref(p);
%
% See also SUBSREF, SUBSASGN, MERGEPREFS

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Developer note:
%
% The reason this object was introduced is to allow for a simplified approach
% to preferences in the upper layers (i.e., CHEBFUN, FUN, and SINGFUN) which
% groups them all into the same namespace while keeping the namespace for
% "tech" preferences separate.  This design mostly achieves this goal with the
% following limited exceptions:
%
%  - Some information really does need to propagate from CHEBFUN all the way
%    down into the tech layer, i.e., CHEBFUN needs to be able to set certain
%    preferences that affect the constructors for the individual techs.
%    Designers of techs should ensure that their classes respond to the
%    following "abstract" preferences in an appropriate manner:  eps,
%    maxLength, exactLength, extrapolate, and sampleTest.
%
%  - The original idea was that the techPrefs field of the CHEBPREF would be
%    the only thing that gets passed to the tech constructor.  This is
%    attainable if one is constructing a CHEBFUN by calling the constructor
%    directly but not, e.g., if one is constructing using COMPOSE().  The
%    reason is that when FUN calls COMPOSE(), it does not know if it is calling
%    the tech COMPOSE() (and therefore should drop the unneeded prefs) or the
%    SINGFUN COMPOSE() (for which it needs to keep them).  Techs can deal with
%    this by dropping the unnecessary information themselves either by doing
%    so directly or by calling CHEBPREF.MERGEPREFS().
%
%    We considered other designs that deal with this problem more gracefully
%    but ultimately found them cumbersome compared to the present one.
%
% The other major design goal was to make this object behave just like a plain
% struct but a little more "intelligently" with respect to system preferences.
% The end result is that CHEBPREFs and ordinary structs are almost
% interchangeable, and one can easily write functions that accept both types of
% arguments by calling P = CHEBPREF(P), where P is the preference input to the
% function.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    properties
        maxTotalLength
        enableBreakpointDetection
        breakpointPrefs
        domain
        enableSingularityDetection
        singPrefs
        tech
        techPrefs
        cheb2Prefs
        % This is a MATLAB structure which stores the system preferences.  As a
        % class invariant, this structure is guaranteed to contain fields for
        % each of the preferences of the upper layers listed above.  It is also
        % guaranteed to contain a techPrefs field, but no guarantees are made
        % about its contents.
        prefList
    end

    methods

        function outPref = chebpref(inPref)
            if ( (nargin == 1) && isa(inPref, 'chebpref') )
                outPref = inPref;
                return
            elseif ( nargin < 1 )
                inPref = struct();
            end

            % Initialize default preference values.
            outPref.maxTotalLength = 65537;
            outPref.enableBreakpointDetection = false;
                outPref.breakpointPrefs.splitMaxLength = 129;
                outPref.breakpointPrefs.splitMaxTotalLength = 6000;
            outPref.domain = [-1 1];
            outPref.enableSingularityDetection = false;
                outPref.singPrefs.exponentTol = 1.1*1e-11;
                outPref.singPrefs.maxPoleOrder = 20;
            outPref.tech = 'chebtech';
            outPref.techPrefs = struct();
                outPref.techPrefs.eps = 2^(-52);
                outPref.techPrefs.maxLength = 65537;
                outPref.techPrefs.exactLength = NaN;
                outPref.techPrefs.extrapolate = false;
                outPref.techPrefs.sampleTest = true;
            outPref.cheb2Prefs = struct(); 
                outPref.cheb2Prefs.maxRank = 2049; 
                outPref.cheb2Prefs.maxLength = 65537; 
                outPref.cheb2Prefs.eps = 2^(-52); 
                outPref.cheb2Prefs.exactLength = false;
                outPref.cheb2Prefs.sampleTest = true; 
            outPref.prefList = chebpref.manageDefaultPrefs('get');

            % Copy fields from q, placing unknown ones in techPrefs and merging
            % incomplete substructures.
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
                    outPref.prefList.techPrefs.(field{1}) = inPref.(field{1});
                end
            end
        end

        function out = subsref(pref, ind)
        %SUBSREF   Subscripted referencing for CHEBPREF.
        %   P.PROP, where P is a CHEBPREF object, returns the value of the
        %   CHEBPREF property PROP stored in P.  If PROP is not a CHEBPREF
        %   property, P.TECHPREFS.PROP will be returned instead.  If PROP is
        %   neither a CHEBPREF property nor a field in P.TECHPREFS, an error
        %   will be thrown.
        %
        %   For access to fields PROP of TECHPREFS that have the same name as a
        %   CHEBPREF property, use the syntax P.TECHPREFS.PROP.
        %
        %   CHEBPREF does not support any other subscripted referencing types,
        %   including '()' and '{}'.
            switch ( ind(1).type )
                case '.'
                    if ( isfield(pref.prefList, ind(1).subs) )
                        out = pref.prefList.(ind(1).subs);
                    else
                        out = pref.prefList.techPrefs.(ind(1).subs);
                    end

                    if ( numel(ind) > 1 )
                        out = subsref(out, ind(2:end));
                    end
                otherwise
                    error('CHEBTECH:subsref:badType', ...
                        'Invalid subscripted reference type.')
            end
        end

        function pref = subsasgn(pref, ind, val)
        %SUBSASGN   Subscripted assignment for CHEBPREF.
        %   P.PROP = VAL, where P is a CHEBPREF object, assigns the value VAL
        %   to the CHEBPREF property PROP stored in P.  If PROP is not a
        %   CHEBPREF property, the assignment will be made to P.TECHPREFS.PROP
        %   instead.
        %
        %   To assign to fields PROP of TECHPREFS that have the same name as a
        %   CHEBPREF property, use the syntax P.TECHPREFS.PROP = VAL.
        %
        %   CHEBPREF does not support any other subscripted assignment types,
        %   including '()' and '{}'.
            switch ( ind(1).type )
                case '.'
                    if ( isfield(pref.prefList, ind(1).subs) )
                        pref.prefList = builtin('subsasgn', pref.prefList, ...
                            ind, val);
                    else
                        pref.prefList.techPrefs = builtin('subsasgn', ...
                            pref.prefList.techPrefs, ind, val);
                    end
                otherwise
                    error('CHEBTECH:subsasgn:badType', ...
                        'Invalid subscripted assignment type.')
            end
        end

        function display(pref)
        %DISPLAY   Display a CHEBPREF object.
        %   DISPLAY(PREF) prints out a list of the preferences stored in the
        %   CHEBPREF object PREF.

            % Compute the screen column in which pref values start.
            valueCol = 34; % length('    enableSingularityDetection:   ');
            for field = fieldnames(pref.prefList.techPrefs).'
                field = field{1};
                col = length(['        ' field '  ']);
                if ( col > valueCol )
                    valueCol = col;
                end
            end

            % A subfunction to pad strings for formatting.
            function s = padString(s)
            %PADSTRING   Add whitespace to string for formatting.
                s = [s repmat(' ', 1, valueCol - length(s))];
            end

            % Print values of "known" preferences.
            prefList = pref.prefList;

            fprintf('chebpref object with the following preferences:\n');
            fprintf([padString('    domain:') '[%g, %g]\n'], ...
                prefList.domain(1), prefList.domain(end));
            fprintf([padString('    enableBreakpointDetection:') '%d\n'], ...
                prefList.enableBreakpointDetection);
            fprintf('    breakpointPrefs\n');
            fprintf([padString('        splitMaxLength:') '%d\n'], ...
                prefList.breakpointPrefs.splitMaxLength');
            fprintf([padString('        splitMaxTotalLength:') '%d\n'], ...
                prefList.breakpointPrefs.splitMaxTotalLength');
            fprintf([padString('    enableSingularityDetection:') '%d\n'], ...
                prefList.enableSingularityDetection);
            fprintf('    singPrefs\n');
            fprintf([padString('        exponentTol:') '%d\n'], ...
                prefList.singPrefs.exponentTol');
            fprintf([padString('        maxPoleOrder:') '%d\n'], ...
                prefList.singPrefs.maxPoleOrder');
            fprintf('    cheb2Prefs\n');
            fprintf([padString('        maxRank:') '%d\n'], ...
                prefList.cheb2Prefs.maxRank');
            fprintf([padString('        maxLength:') '%d\n'], ...
                prefList.cheb2Prefs.maxLength');            
            fprintf([padString('        eps:') '%d\n'], ...
                prefList.cheb2Prefs.eps');            
            fprintf([padString('        exactLength:') '%d\n'], ...
                prefList.cheb2Prefs.exactLength');            
            fprintf([padString('        sampleTest:') '%d\n'], ...
                prefList.cheb2Prefs.sampleTest');            
            fprintf([padString('    scale:') '%d\n'], ...
                prefList.scale);
            fprintf([padString('    tech:') '''%s''\n'], ...
                prefList.tech)
            fprintf('    techPrefs\n');

            % Format and print values of tech preferences.
            for field = fieldnames(prefList.techPrefs).'
                field = field{1};
                printStr = padString(['        ' field ':']);

                if ( isempty(prefList.techPrefs.(field)) )
                    fprintf([printStr 'empty\n']);
                elseif ( ischar(prefList.techPrefs.(field)) && ...
                         isrow(prefList.techPrefs.(field)) )
                    fprintf([printStr '''%s''\n'], prefList.techPrefs.(field))
                elseif ( numel(prefList.techPrefs.(field)) > 1 )
                    fprintf([printStr class(prefList.techPrefs.(field)) ...
                        ' array\n']);
                elseif ( isfloat(prefList.techPrefs.(field)) )
                    fprintf([printStr '%0.16g\n'], prefList.techPrefs.(field))
                elseif ( islogical(prefList.techPrefs.(field)) )
                    fprintf([printStr '%d\n'], prefList.techPrefs.(field))
                else
                    fprintf([printStr class(prefList.techPrefs.(field)) '\n']);
                end
            end
        end

    end

    methods ( Static = true )

        function pref1 = mergePrefs(pref1, pref2, map)
        %MERGEPREFS   Merge preference structures.
        %   P = CHEBPREF.MERGEPREFS(P, Q), where P and Q are MATLAB structures,
        %   "merges" Q into P by replacing the contents of fields in P with
        %   those of identically-named fields in Q.  If Q has a field whose
        %   name does not match any of those in P, it is added to P.
        %
        %   P = CHEBPREF.MERGEPREFS(P, Q, MAP) does the same but uses the
        %   structure MAP to "translate" the names of fields of Q into names of
        %   fields of P.  If Q has a field FIELD and MAP has a field of the
        %   same name, then the value of P.(MAP.FIELD) will be replaced by the
        %   contents of Q.FIELD.  If P does not have a field matching the
        %   string stored in MAP.FIELD, one will be added to P.
        %
        %   P and Q may also be CHEBPREF objects.  In this case, P and Q are
        %   replaced by P.TECHPREFS and Q.TECHPREFS before proceeding, and the
        %   output is a MATLAB structure suitable for passing as a preference
        %   argument to a "tech" constructor.

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Developer notes:
        %  - This function is a helper function intended for use by "technology"
        %    objects (usually subclasses of SMOOTHFUN) for managing their
        %    preferences.  See CHEBTECH.TECHPREF for an illustration.
        %  - The second syntax is useful, e.g., if Q contains abstractly-named
        %    preferences which may have a better name within the specific
        %    context of the tech object whose preferences are stored in P.
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            if ( isa(pref1, 'chebpref') )
                pref1 = pref1.prefList.techPrefs;
            end

            if ( isa(pref2, 'chebpref') )
                pref2 = pref2.prefList.techPrefs;
            end

            if ( nargin < 3 )
                map = struct();
            end

            for field = fieldnames(pref2).'
                if ( isfield(map, field{1}) )
                    pref1.(map.(field{1})) = pref2.(field{1});
                else
                    pref1.(field{1}) = pref2.(field{1});
                end
            end
        end

        function pref = getFactoryDefaults(getFactory)
        %GETFACTORYDEFAULTS   Get factory default preferences.
        %   PREF = CHEBPREF.GETFACTORYDEFAULTS() returns a CHEBPREF object with
        %   the preferences set to their factory defaults, irrespective of the
        %   currently defined values of the default preferences.  This function
        %   is useful if the user wishes to construct a CHEBFUN using the
        %   factory defaults when other user-set defaults are currently in
        %   force.
        %
        % See also GETDEFAULTS, SETDEFAULTS.

            fd = chebpref.factoryDefaultPrefs();
            pref = chebpref(fd);

            % Undo any merging of the factory default techPrefs with the
            % currently defined default techPrefs.  This is necessary, e.g., if
            % the current defaults have techPrefs stored that are not among the
            % factory defaults.
            pref.prefList.techPrefs = fd.techPrefs;
        end

        function pref = getDefaults()
        %GETDEFAULTS   Get default preferences.
        %   PREF = CHEBPREF.GETDEFAULTS() returns a CHEBPREF object with the
        %   preferences set to the currently stored default values.  It is
        %   equivalent to PREF = CHEBPREF().
        %
        % See also GETFACTORYDEFAULTS, SETDEFAULTS.

            pref = chebpref();
        end

        function setDefaults(varargin)
        %SETDEFAULTS   Set default preferences.
        %   CHEBPREF.SETDEFAULTS(PREF1, VAL1, PREF2, VAL2, ...) sets the
        %   default values for the preferences whose names are stored in the
        %   strings PREF1, PREF2, ..., etc. to VAL1, VAL2, ..., etc.  All
        %   subsequently constructed CHEBPREF objects will use these values as
        %   the defaults.
        %
        %   CHEBPREF.SETDEFAULTS(PREF) sets the default values to the
        %   preferences stored in the CHEBPREF object PREF.  PREF can also be a
        %   MATLAB structure, in which case it is converted to a CHEBPREF as
        %   described in the documentation for the CHEBPREF constructor first.
        %
        %   CHEBPREF.SETDEFAULTS('factory') resets the default preferences to
        %   their factory values.
        %
        % See also GETDEFAULTS, GETFACTORYDEFAULTS.

        % TODO:  What to do about preferences stored in substructures, like
        % singfun.exponentTol?  Aside from preferences in techPrefs whose names
        % don't collide with other "top-level" preferences, these can't be set
        % using the first syntax listed above (though they still can be set
        % with the second).

            if ( nargin < 1)
                error('CHEBPREF:setDefaults:notEnoughArguments', ...
                    'Not enough arguments.');
            end

            if ( nargin == 1 )
                if ( isstruct(varargin{1}) )
                    varargin{1} = chebpref(varargin{1});
                end

                if ( ischar(varargin{1}) && strcmp(varargin{1}, 'factory') )
                    chebpref.manageDefaultPrefs('set-factory');
                elseif ( isa(varargin{1}, 'chebpref') )
                    chebpref.manageDefaultPrefs('set', varargin{1}.prefList);
                else
                    error('CHEBPREF:setDefaults:badArg', ...
                        ['When calling chebpref.setDefaults() with just ' ...
                         'one argument, that argument must be ''factory'', ' ...
                         'a CHEBPREF object, or a MATLAB structure.']);
                end
            elseif ( mod(nargin, 2) == 0 )
                chebpref.manageDefaultPrefs('set', varargin{:});
            else
                error('CHEBPREF:setDefaults:unpairedArg', ...
                    'Unpaired argument in name-value pair list.');
            end
        end
    end

    methods ( Static = true, Access = private)

        function varargout = manageDefaultPrefs(varargin)
        %MANAGEDEFAULTPREFS   Private method for handling default preferences.
        %   CHEBPREF.MANAGEDEFAULTPREFS('get') returns a structure suitable for
        %   storing in the prefList property of a CHEBPREF with all of the
        %   currently stored default preferences suitable for initializing a
        %   CHEBPREF object.
        %
        %   CHEBPREF.MANAGEDEFAULTPREFS('set-factory') restores the default
        %   preferences to their "factory" values.
        %
        %   CHEBPREF.MANAGEDEFAULTPREFS('set', PREFLIST) sets the default
        %   values to those stored in the structure PREFLIST.  PREFLIST should
        %   be a structure suitable for use as a CHEBPREF prefList.
        %
        %   CHEBPREF.MANAGEDEFAULTPREFS('set', PREF1, VAL1, PREF2, VAL2, ...)
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
        %    CHEBPREF as a MATLAB structure, and that's not something with which
        %    anyone outside of CHEBPREF should be concerned.
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            persistent defaultPrefs;

            if ( isempty(defaultPrefs) )
                defaultPrefs = chebpref.factoryDefaultPrefs();
            end

            if ( strcmp(varargin{1}, 'get') )
                varargout{1} = defaultPrefs;
            elseif ( strcmp(varargin{1}, 'set-factory') )
                defaultPrefs = chebpref.factoryDefaultPrefs();
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
                            defaultPrefs.techPrefs.(prefName) = prefValue;
                        end
                        varargin(1:2) = [];
                    end
                end
            end
        end

        function factoryPrefs = factoryDefaultPrefs()
        %FACTORYDEFAULTPREFS   Get structure of factory default preferences.
        %   S = CHEBPREF.FACTORYDEFAULTREFS() returns a structure suitable for
        %   storing in the prefList property of a CHEBPREF object that contains
        %   all of the "factory default" values of the CHEBFUN
        %   construction-time preferences.

            factoryPrefs.domain = [-1 1];
            factoryPrefs.enableBreakpointDetection = false;
                factoryPrefs.breakpointPrefs.splitMaxLength = 129;
                factoryPrefs.breakpointPrefs.splitMaxTotalLength = 6000;
            factoryPrefs.enableSingularityDetection = false;
                factoryPrefs.singPrefs.exponentTol = 1.1*1e-11;
                factoryPrefs.singPrefs.maxPoleOrder = 20;
                factoryPrefs.singPrefs.exponents = [];
                factoryPrefs.singPrefs.singType = {};
            factoryPrefs.enableDeltaFunctions = false;
                factoryPrefs.deltaPrefs.deltaTol = 1e-11;
                factoryPrefs.deltaPrefs.proximityTol = 1e-11;

            factoryPrefs.scale = 0;
            factoryPrefs.tech = 'chebtech';
            factoryPrefs.techPrefs = struct();
                factoryPrefs.techPrefs.eps = 2^(-52);
                factoryPrefs.techPrefs.maxLength = 65537;
                factoryPrefs.techPrefs.exactLength = NaN;
                factoryPrefs.techPrefs.extrapolate = false;
                factoryPrefs.techPrefs.sampleTest = true;
            factoryPrefs.cheb2Prefs = struct(); 
                factoryPrefs.cheb2Prefs.maxRank = 2049;   
                factoryPrefs.cheb2Prefs.maxLength = 65537;   
                factoryPrefs.cheb2Prefs.eps = 2^(-52);   
                factoryPrefs.cheb2Prefs.exactLength = 0; 
                factoryPrefs.cheb2Prefs.sampleTest = 1;
        end

    end
end

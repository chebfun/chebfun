classdef chebfunpref < chebpref
%CHEBFUNPREF   Class for managing Chebfun construction-time preferences.
%   CHEBFUNPREF is a class for managing Chebfun construction-time preferences
%   such as the construction tolerance, whether or not to perform breakpoint and
%   singularity detection, and the various options that those features require.
%   These objects can be supplied to the CHEBFUN constructor (as well as the
%   constructors of other classes in Chebfun), which will interpret them and
%   adjust the construction process accordingly. Note that all preferences _are_
%   case sensitive.
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
%   splitting                  - Enable/disable breakpoint detection.
%     true
%    [false]
%
%     If true, breakpoints between FUNS may be introduced where a discontinuity
%     in a function or a low-order derivative is detected or if a global
%     representation will be too long.  If false, breakpoints will be
%     introduced only at points where discontinuities are being created (e.g.,
%     by ABS(F) at points where a CHEBFUN F passes through zero).
%
%   splitPrefs                 - Preferences for breakpoint detection.
%
%      splitLength             - Maximum FUN length.
%       [160]
%
%         This is the maximum length of a single FUN (e.g., the number of
%         Chebyshev points used for FUNs based on Chebyshev polynomial
%         interpolation) allowed by the constructor when breakpoint detection
%         is enabled.
%
%      splitMaxLength          - Maximum total CHEBFUN length.
%       [6000]
%
%         This is the maximum total length of the CHEBFUN (i.e., the sum of the
%         lengths of all the FUNs) allowed by the constructor when breakpoint
%         detection is enabled.
%
%   blowup                     - Enable/disable singularity detection.
%     true
%    [false]
%
%      If true, the constructor will attempt to detect and factor out
%      singularities, (e.g., points where a function or its derivatives become
%      unbounded). If false, breakpoints will be introduced only at points where
%      singularities are being created, (e.g., by SQRT(F) at points where a
%      CHEBFUN F passes through zero). See SINGFUN for more information.
%
%   blowupPrefs                - Preferences for blowup / singularity detection.
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
%      defaultSingType         - Type of singularities.
%         
%         The default singularity type to be used when singularity detection is
%         enabled and no singType is provided.
%
%   enableDeltaFunctions - Enable delta functions.
%     true
%    [false]
%       
%      If true, the DELTAFUN class will be invoked to manage any delta 
%      functions present in the object.
%
%   deltaPrefs                 - Preferences for delta functions.
%
%      deltaTol                - Tolerance for magnitude of delta functions.
%       [1e-11]
%
%         This is the tolerance up to which delta functions will be negligible.
%
%      proximityTol            - Minimum distance between delta functions.
%       [1e-11]
%
%         If two delta functions are located closer than this tolerance, they 
%         will be merged.
%
%   tech                       - Representation technology.
%    ['chebtech2']
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
%      fixedLength             - Exact representation length.
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
% using CHEBFUNPREF.SETDEFAULTS(); see the documentation for that function for
% further details.
%
% Constructor inputs:
%   P = CHEBFUNPREF() creates a CHEBFUNPREF object with the default values of
%   the preferences.  For a list of all available preferences, see above.
%
%   P = CHEBFUNPREF(Q), where Q is a MATLAB structure uses the field/value pairs
%   in Q to set the properties of P.  If a field of Q has a name which matches
%   a property of P, the value of that property of P is set to the value
%   associated to that field in Q.  Any fields of Q that are not properties of
%   P are interpreted as preferences for the constructor of the underlying
%   representation technology and are placed in P.TECHPREFS.  The exceptions to
%   this are the fields SPLITPREFS, BLOWUPPREFS, and TECHPREFS.  If Q has
%   fields with these names, they will be assumed to be MATLAB structures and
%   will be "merged" with the structures of default preferences stored in the
%   properties of the same names in P using CHEBFUNPREF.MERGEPREFS().
%
%   P = CHEBFUNPREF(Q), where Q is a CHEBFUNPREF, sets P to be a copy of Q.
%
%   R = CHEBFUNPREF(P, Q), where P is a CHEBFUNPREF and Q is a MATLAB
%   structure, is similar to CHEBFUNPREF(Q) except that the preferences in P
%   are used as the base set of preferences instead of the currently stored
%   defaults.  The output R will be a CHEBFUNPREF with the preferences of P
%   overridden by the field/value pairs in the structure Q.
%
%   R = CHEBFUNPREF(P, Q), where P and Q are both CHEBFUNPREF objects is
%   similar to the previous syntax.  The output R is a CHEBFUNPREF with the
%   preferences of P overridden by those in Q.  This is equivalent to setting R
%   to be a copy of Q plus any additional TECHPREFS stored in P that were not
%   stored in Q.
%
% Notes: Creating preferences from structures.
%   When building a CHEBFUNPREF from a structure using the second calling
%   syntax above, one should take care to ensure that preferences for the
%   underlying representation technology are specified once and only once;
%   e.g., do not simultaneously set Q.MYPREF = 1 and Q.TECHPREFS.MYPREF = 2.
%   The value of P.TECHPREFS.MYPREF that gets set from P = CHEBFUNPREF(Q) in
%   this circumstance is implementation-defined.
%
% Notes: Relationship to V4 preferences.
%   All V4 preference names are supported in calls to the CHEBFUN constructor
%   (although support for this may be removed in a future release). Here is a
%   rough guide to how the new V5 preferences relate to the old V4 ones. Note
%   that the V5 names _must_ be used in CHEBFUNPREF(), which does _not_ support
%   the V4 names.
%        V4                  V5
%    'maxdegree'   -->   'maxLength'
%    'maxlength'   -->  {'splitPrefs', 'splitMaxLength'}
%    'splitdegree' -->  {'splitPrefs', 'splitLength'}
%    'resampling'  -->  'refinementFunction'
%    'sampletest'  -->  'sampleTest'
%    'chebkind'    -->  'tech'
%    'plot_numpts'       removed
%    'polishroots'       removed
%    'ADdepth'           removed
%   (Note that when setting preferences directly via the constructor, one should
%   only include the second entry in those preferences given as cell arrays
%   above. For example, chebfun(@abs, 'splitLength', 10);.)
%
% Examples:
%   Create a CHEBFUNPREF for building a CHEBFUN based on CHEBTECH (default) with
%   breakpoint detection, a splitting length of 257 (pieces of polynomial degree
%   256, and a custom CHEBTECH refinement function:
%      p.splitting = true;
%      p.splitPrefs.splitLength = 257;
%      p.techPrefs.refinementFunction = @custom;
%      pref = chebfunpref(p);
%
%   Same thing with a slightly shorter syntax:
%      p.splitting = true;
%      p.splitPrefs.splitLength = 257;
%      p.refinementFunction = @custom;
%      pref = chebfunpref(p);
%
% See also CHEBOPPREF.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEVELOPER NOTE:
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
%    maxLength, fixedLength, extrapolate, and sampleTest.
%
%  - The original idea was that the techPrefs field of the CHEBFUNPREF would be
%    the only thing that gets passed to the tech constructor.  This is
%    attainable if one is constructing a CHEBFUN by calling the constructor
%    directly but not, e.g., if one is constructing using COMPOSE().  The
%    reason is that when FUN calls COMPOSE(), it does not know if it is calling
%    the tech COMPOSE() (and therefore should drop the unneeded prefs) or the
%    SINGFUN COMPOSE() (for which it needs to keep them).  Techs can deal with
%    this by dropping the unnecessary information themselves either by doing
%    so directly or by calling CHEBFUNPREF.MERGEPREFS().
%
%    We considered other designs that deal with this problem more gracefully
%    but ultimately found them cumbersome compared to the present one.
%
% The other major design goal was to make this object behave just like a plain
% struct but a little more "intelligently" with respect to system preferences.
% The end result is that CHEBFUNPREFs and ordinary structs are almost
% interchangeable, and one can easily write functions that accept both types of
% arguments by calling P = CHEBFUNPREF(P), where P is the preference input to
% the function.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CLASS CONSTRUCTOR:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods ( Access = public, Static = false )

        function outPref = chebfunpref(varargin)
            if ( nargin < 1 )
                inPrefList = struct();
            elseif ( ischar(varargin{1}) )
                if ( nargin == 1 )
                    error('CHEBFUN:CHEBFUNPREF:chebfunpref:deprecated', ...
                        ['chebfunpref() no longer supports queries of ', ...
                         'the form chebfunpref(''prop'').\n', ...
                         'Please use chebfunpref().prop.']);
                else
                    error('CHEBFUN:CHEBFUNPREF:chebfunpref:deprecated', ...
                        ['chebfunpref() no longer supports assignment ', ...
                         'via chebfunpref(''prop'', val).\n', ...
                         'Please use chebfunpref.setDefaults(''prop'', val).']);
                end
            elseif ( nargin == 1 )
                if ( isa(varargin{1}, 'chebfunpref') )
                    outPref = varargin{1};
                    return
                else
                    inPrefList = varargin{1};
                end
            elseif ( nargin == 2 )
                if ( ~isa(varargin{1}, 'chebfunpref') || ...
                     (~isa(varargin{2}, 'chebfunpref') && ...
                      ~isstruct(varargin{2})) )
                      error('CHEBFUN:CHEBFUNPREF:chebfunpref:badTwoArgCall', ...
                        ['When calling CHEBFUNPREF with two arguments, ' ...
                        'the first must be a CHEBFUNPREF, and the ' ...
                        'second must be a CHEBFUNPREF or a struct.']);
                elseif ( isa(varargin{2}, 'chebfunpref') )
                    inPrefList = varargin{2}.prefList;
                else
                    inPrefList = varargin{2};
                end
            elseif ( nargin > 2 )
                error('CHEBFUN:CHEBFUNPREF:chebfunpref:tooManyInputs', ...
                    'Too many input arguments.')
            end

            % If the user supplied a base set of preferences to be overridden,
            % use it; otherwise, use the stored defaults.
            if ( nargin <= 1 )
                outPref.prefList = chebfunpref.manageDefaultPrefs('get');
            else
                outPref.prefList = varargin{1}.prefList;
            end

            % Copy fields from inPrefList, placing unknown ones in techPrefs
            % and merging incomplete substructures.
            for field = fieldnames(inPrefList).'
                if ( isfield(outPref.prefList, field{1}) )
                    if ( isstruct(outPref.prefList.(field{1})) )
                        outPref.prefList.(field{1}) = ...
                            chebfunpref.mergeTechPrefs(...
                                outPref.prefList.(field{1}), ...
                                inPrefList.(field{1}));
                    else
                        outPref.prefList.(field{1}) = inPrefList.(field{1});
                    end
                else
                    outPref.prefList.techPrefs.(field{1}) = ...
                        inPrefList.(field{1});
                end
            end
        end
        
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CLASS METHODS:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods ( Access = public, Static = false )

        function out = subsref(pref, ind)
        %SUBSREF   Subscripted referencing for CHEBFUNPREF.
        %   P.PROP, where P is a CHEBFUNPREF object, returns the value of the
        %   CHEBFUNPREF property PROP stored in P.  If PROP is not a CHEBFUNPREF
        %   property, P.TECHPREFS.PROP will be returned instead.  If PROP is
        %   neither a CHEBFUNPREF property nor a field in P.TECHPREFS, an error
        %   will be thrown.
        %
        %   For access to fields PROP of TECHPREFS that have the same name as a
        %   CHEBFUNPREF property, use the syntax P.TECHPREFS.PROP.
        %
        %   CHEBFUNPREF does not support any other subscripted referencing
        %   types, including '()' and '{}'.

            switch ( ind(1).type )
                case '.'
                    if ( isfield(pref.prefList, ind(1).subs) )
                        out = pref.prefList.(ind(1).subs);
                    else
                        techObj = feval(pref.prefList.tech);
                        fullTechPrefs = ...
                            techObj.techPref(pref.prefList.techPrefs);
                        if ( isfield(fullTechPrefs, ind(1).subs) )
                            % Try to find the tech preference name after
                            % merginging with the full list of tech preferences
                            % obtained via the tech's techPref() function of the
                            % current tech.
                            out = fullTechPrefs.(ind(1).subs);
                        else
                            % If we couldn't find the tech preference name
                            % above, it may be because it was an abstractly
                            % named preference that got mapped to something the
                            % tech's techPref() deemed more sensible.  So, we
                            % also try looking in the list of tech preferences
                            % we have prior to forming the full list.
                            out = pref.prefList.techPrefs.(ind(1).subs);
                        end
                    end

                    if ( numel(ind) > 1 )
                        out = subsref(out, ind(2:end));
                    end
                otherwise
                    error('CHEBFUN:CHEBFUNPREF:subsref:badType', ...
                        'Invalid subscripted reference type.')
            end
        end

        function pref = subsasgn(pref, ind, val)
        %SUBSASGN   Subscripted assignment for CHEBFUNPREF.
        %   P.PROP = VAL, where P is a CHEBFUNPREF object, assigns the value
        %   VAL to the CHEBFUNPREF property PROP stored in P.  If PROP is not a
        %   CHEBFUNPREF property, the assignment will be made to
        %   P.TECHPREFS.PROP instead.
        %
        %   To assign to fields PROP of TECHPREFS that have the same name as a
        %   CHEBFUNPREF property, use the syntax P.TECHPREFS.PROP = VAL.
        %
        %   CHEBFUNPREF does not support any other subscripted assignment types,
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
                    error('CHEBFUN:CHEBFUNPREF:subsasgn:badType', ...
                        'Invalid subscripted assignment type.')
            end
        end

        function display(pref)
        %DISPLAY   Display a CHEBFUNPREF object.
        %   DISPLAY(PREF) prints out a list of the preferences stored in the
        %   CHEBFUNPREF object PREF.

            prefList = pref.prefList;

            % Get a list of the merged tech preferences.
            tech = prefList.tech;
            techObj = feval(tech);
            techPrefs = techObj.techPref(prefList.techPrefs);

            % Compute the screen column in which pref values start.
            valueCol = 34; % length('    blowup:   ');
            for field = fieldnames(techPrefs).'
                field1 = field{1};
                col = length(['        ' field1 '  ']);
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
            fprintf('chebfunpref object with the following preferences:\n');
            fprintf([padString('    domain:') '[%g, %g]\n'], ...
                prefList.domain(1), prefList.domain(end));
            fprintf([padString('    splitting:') '%d\n'], ...
                prefList.splitting);
            fprintf('    splitPrefs\n');
            fprintf([padString('        splitLength:') '%d\n'], ...
                prefList.splitPrefs.splitLength');
            fprintf([padString('        splitMaxLength:') '%d\n'], ...
                prefList.splitPrefs.splitMaxLength');
            fprintf([padString('    blowup:') '%d\n'], ...
                prefList.blowup);
            fprintf('    blowupPrefs\n');
            fprintf([padString('        exponentTol:') '%d\n'], ...
                prefList.blowupPrefs.exponentTol');
            fprintf([padString('        maxPoleOrder:') '%d\n'], ...
                prefList.blowupPrefs.maxPoleOrder');
            fprintf([padString('        defaultSingType:') '''%s''\n'], ...
                prefList.blowupPrefs.defaultSingType');            
            fprintf([padString('    enableDeltaFunctions:') '%d\n'], ...
                prefList.enableDeltaFunctions);
            fprintf('    deltaPrefs\n');
            fprintf([padString('        deltaTol:') '%d\n'], ...
                prefList.deltaPrefs.deltaTol');
            fprintf([padString('        proximityTol:') '%d\n'], ...
                prefList.deltaPrefs.proximityTol');                      
            fprintf('    cheb2Prefs\n');
            fprintf([padString('        maxRank:') '%d\n'], ...
                prefList.cheb2Prefs.maxRank');
            fprintf([padString('        sampleTest:') '%d\n'], ...
                prefList.cheb2Prefs.sampleTest');
            
            techStr = func2str(tech);
            fprintf([padString('    tech:') '@%s\n'], techStr)

            if ( usejava('jvm') && usejava('desktop') )
                link = '    <a href="matlab: help %s/techPref">techPrefs</a>\n';
                fprintf(link, techStr)
            else
                fprintf('    techPrefs\n')
            end

            % Format and print values of tech preferences.
            for field = fieldnames(techPrefs).'
                field1 = field{1};
                printStr = padString(['        ' field1 ':']);

                if ( isempty(techPrefs.(field1)) )
                    fprintf([printStr 'empty\n']);
                elseif ( ischar(techPrefs.(field1)) && ...
                         isrow(techPrefs.(field1)) )
                    fprintf([printStr '''%s''\n'], techPrefs.(field1))
                elseif ( numel(techPrefs.(field1)) > 1 )
                    fprintf([printStr class(techPrefs.(field1)) ...
                        ' array\n']);
                elseif ( isfloat(techPrefs.(field1)) )
                    fprintf([printStr '%0.16g\n'], techPrefs.(field1))
                elseif ( islogical(techPrefs.(field1)) )
                    fprintf([printStr '%d\n'], techPrefs.(field1))
                else
                    fprintf([printStr class(techPrefs.(field1)) '\n']);
                end
            end
        end

    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% STATIC METHODS:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    methods ( Access = public, Static = true )

        function pref1 = mergeTechPrefs(pref1, pref2)
        %MERGETECHPREFS   Merge tech preference structures.
        %   P = CHEBFUNPREF.MERGETECHPREFS(P, Q), where P and Q are MATLAB
        %   structures, "merges" Q into P by replacing the contents of fields
        %   in P with those of identically-named fields in Q.  If Q has a field
        %   whose name does not match any of those in P, it is added to P.
        %
        %   P and Q may also be CHEBFUNPREF objects.  In this case, P and Q are
        %   replaced by P.TECHPREFS and Q.TECHPREFS before proceeding, and the
        %   output is a MATLAB structure suitable for passing as a preference
        %   argument to a "tech" constructor.

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % DEVELOPER NOTES:
        %  - This function is a helper function intended for use by "technology"
        %    objects (usually subclasses of SMOOTHFUN) for managing their
        %    preferences.  See CHEBTECH.TECHPREF for an illustration.
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            if ( isa(pref1, 'chebfunpref') )
                pref1 = pref1.prefList.techPrefs;
            end

            if ( isa(pref2, 'chebfunpref') )
                pref2 = pref2.prefList.techPrefs;
            end

            pref1 = chebpref.mergePrefStructs(pref1, pref2);
        end

        function pref = getFactoryDefaults()
        %GETFACTORYDEFAULTS   Get factory default preferences.
        %   PREF = CHEBFUNPREF.GETFACTORYDEFAULTS() returns a CHEBFUNPREF
        %   object with the preferences set to their factory defaults,
        %   irrespective of the currently defined values of the default
        %   preferences.  This function is useful if the user wishes to
        %   construct a CHEBFUN using the factory defaults when other user-set
        %   defaults are currently in force.
        %
        % See also SETDEFAULTS.

            fd = chebfunpref.factoryDefaultPrefs();
            pref = chebfunpref(fd);

            % Undo any merging of the factory default techPrefs with the
            % currently defined default techPrefs.  This is necessary, e.g., if
            % the current defaults have techPrefs stored that are not among the
            % factory defaults.
            pref.prefList.techPrefs = fd.techPrefs;
        end

        function setDefaults(varargin)
        %SETDEFAULTS   Set default preferences.
        %   CHEBFUNPREF.SETDEFAULTS(PREF1, VAL1, PREF2, VAL2, ...) sets the
        %   default values for the preferences whose names are stored in the
        %   strings PREF1, PREF2, ..., etc. to VAL1, VAL2, ..., etc.  All
        %   subsequently constructed CHEBFUNPREF objects will use these values
        %   as the defaults.
        %
        %   To set defaults for second tier preferences, such as
        %   splitPrefs.splitLength, one can use the syntax
        %   CHEBFUNPREF.SETDEFAULT({'splitPrefs', 'splitLength'}, 257).
        %   However, this syntax is still experimental.
        %
        %   CHEBFUNPREF.SETDEFAULTS(PREF) sets the default values to the
        %   preferences stored in the CHEBFUNPREF object PREF.  PREF can also
        %   be a MATLAB structure, in which case it is converted to a
        %   CHEBFUNPREF as described in the documentation for the CHEBFUNPREF
        %   constructor first.
        %
        %   CHEBFUNPREF.SETDEFAULTS('factory') resets the default preferences to
        %   their factory values.
        %
        % See also GETFACTORYDEFAULTS.

        % TODO:  What to do about preferences stored in substructures, like
        % singfun.exponentTol?  Aside from preferences in techPrefs whose names
        % don't collide with other "top-level" preferences, these can't be set
        % using the first syntax listed above (though they still can be set
        % with the second).  

            % The reason we don't just use manageDefaults as the second
            % argument to chebpref.setDefaults and wrap it in an additional
            % anonymous function instead is to get around what seems to be a
            % bug in MATLAB.  See commit messages for more information.
            manageDefaults = @chebfunpref.manageDefaultPrefs;
            chebpref.setDefaults(@(inPref) chebfunpref(inPref), ...
                @(varargin) manageDefaults(varargin{:}), varargin{:});
        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% PRIVATE STATIC METHODS
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods ( Static = true, Access = private )

        function varargout = manageDefaultPrefs(varargin)
        %MANAGEDEFAULTPREFS   Private method for handling default preferences.
        %   CHEBFUNPREF.MANAGEDEFAULTPREFS('get') returns a structure suitable
        %   for storing in the prefList property of a CHEBFUNPREF with all of
        %   the currently stored default preferences suitable for initializing
        %   a CHEBFUNPREF object.
        %
        %   CHEBFUNPREF.MANAGEDEFAULTPREFS('set-factory') restores the default
        %   preferences to their "factory" values.
        %
        %   CHEBFUNPREF.MANAGEDEFAULTPREFS('set', PREFLIST) sets the default
        %   values to those stored in the structure PREFLIST.  PREFLIST should
        %   be a structure suitable for use as a CHEBFUNPREF prefList.
        %
        %   CHEBFUNPREF.MANAGEDEFAULTPREFS('set', PREF1, VAL1, PREF2, VAL2, ...)
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
        %    CHEBFUNPREF as a MATLAB structure, and that's not something with
        %    which anyone outside of CHEBFUNPREF should be concerned.
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            persistent defaultPrefs;

            if ( isempty(defaultPrefs) )
                defaultPrefs = chebfunpref.factoryDefaultPrefs();
            end

            if ( strcmp(varargin{1}, 'get') )
                varargout{1} = defaultPrefs;
            elseif ( strcmp(varargin{1}, 'set-factory') )
                defaultPrefs = chebfunpref.factoryDefaultPrefs();
            elseif ( strcmp(varargin{1}, 'set') )
                    varargin(1) = [];
                if ( isstruct(varargin{1}) )
                    defaultPrefs = varargin{1};
                else
                    while ( ~isempty(varargin) )
                        prefName = varargin{1};
                        prefValue = varargin{2};
                        if ( iscell(prefName) && ...
                                isfield(defaultPrefs.(prefName{1}), prefName{2}) )
                            % TODO: Revisit the syntax for assigning to second
                            % tier preferences.
                            defaultPrefs.(prefName{1}).(prefName{2}) = prefValue;
                        elseif ( isfield(defaultPrefs, prefName) )
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
        %   S = CHEBFUNPREF.FACTORYDEFAULTPREFS() returns a structure suitable
        %   for storing in the prefList property of a CHEBFUNPREF object that
        %   contains all of the "factory default" values of the CHEBFUN
        %   construction-time preferences.

            factoryPrefs.domain = [-1 1];
            factoryPrefs.splitting = false;
                factoryPrefs.splitPrefs.splitLength = 160;
                factoryPrefs.splitPrefs.splitMaxLength = 6000;
            factoryPrefs.blowup = false;
                factoryPrefs.blowupPrefs.exponentTol = 1.1*1e-11;
                factoryPrefs.blowupPrefs.maxPoleOrder = 20;
                factoryPrefs.blowupPrefs.defaultSingType = 'sing';                
            factoryPrefs.enableDeltaFunctions = true;
                factoryPrefs.deltaPrefs.deltaTol = 1e-9;
                factoryPrefs.deltaPrefs.proximityTol = 1e-11;
            factoryPrefs.tech = @chebtech2;
            factoryPrefs.techPrefs = struct();
                factoryPrefs.techPrefs.eps = 2^(-52);
                factoryPrefs.techPrefs.maxLength = 65537;
                factoryPrefs.techPrefs.fixedLength = NaN;
                factoryPrefs.techPrefs.extrapolate = false;
                factoryPrefs.techPrefs.sampleTest = true;
            factoryPrefs.cheb2Prefs = struct(); 
                factoryPrefs.cheb2Prefs.maxRank = 513;   
                factoryPrefs.cheb2Prefs.sampleTest = 1;
        end

    end
    
end

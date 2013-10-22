classdef chebpref
%CHEBPREF   Class for managing Chebfun preferences
%   CHEBPREF is a class for managing Chebfun construction-time preferences such
%   as the construction tolerance, whether or not to perform breakpoint and
%   singularity detection, and the various options that those features require.
%   These objects can be supplied to the CHEBFUN constructor, which will
%   interpret them and adjust the construction process accordingly.
%
% Available Preferences:
%
%   maxTotalLength              - Maximum total CHEBFUN length.
%    [65536]
%
%      Sets the maximum allowed "length" of the constructed CHEBFUN when
%      breakpoint detection is disabled, where the notion of length is
%      determined by the underlying representation technology (e.g., polynomial
%      degree, for Chebyshev polynomial interpolation).
%
%   enableBreakpointDetection  - Enable/disable breakpoint detection.
%     true
%    [false]
%
%     If true, breakpoints between funs may be introduced where a discontinuity
%     in a function or a low-order derivative is detected or if a global
%     representation will be too long.  If false, breakpoints will be
%     introduced only at points where discontinuities are being created, e.g.,
%     by ABS(F) at points where a CHEBFUN F passes through zero.
%
%   breakpointPrefs            - Preferences for breakpoint detection.
%
%      splitMaxLength          - Maximum FUN length.
%       [128]
%
%         This is the maximum length of a single FUN (i.e., polynomial degree
%         for funs based on Chebyshev polynomial representations) allowed by
%         the constructor when breakpoint detection is enabled.
%
%      splitMaxTotalLength     - Maximum total CHEBFUN length.
%       [6000]
%
%         This is the maximum total length of the CHEBFUN (i.e., the sum of the
%         lengths of all the FUNs) allowed by the constructor when breakpoint
%         detection is enabled.
%
%   domain                     - Construction domain.
%    [-1, 1]
%
%      This sets the default domain that will be used for CHEBFUN and/or FUN
%      construction if no domain argument is explicitly passed to the
%      constructor.
%
%   enableSingularityDetection - Enable/disable singularity detection.
%     true
%    [false]
%
%      If true, the constructor will attempt to detect and factor out
%      singularities, e.g., points where a function or its derivatives become
%      unbounded.  See SINGFUN for more information.
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
%   enableFunqui               - Enable/disable FUNQUI.
%     true
%    [false]
%
%      Enables the use of the FUNQUI pre-processor for stable interpolation of
%      data from a grid of equally-spaced points.  See documentation for the
%      'equi' flag for CHEBFUN for more information.
%
%   tech                       - Representation technology.
%    ['chebtech']
%
%      Sets the underlying representation technology used to construct the FUNs.
%
%   techPrefs                  - Preferences for the tech constructor.
%
%      This is a structure of preferences that will be passed to the constructor
%      for the underlying representation technology.  See CHEBTECH/PREF for
%      preferences accepted by the default CHEBTECH technology.  Additionally,
%      all techs are required to accept the following preferences:
%
%      eps                     - Construction tolerance.
%       [2^(-52]
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
%        Exact length of the underlying representation.  A NaN value indicates
%        that any length (up to maxLength) is permissible.
%
%      extrapolate             - Extrapolate endpoint values.
%        true
%       [false]
%
%        If true, the tech should avoid direct evaluation of the function at
%        the interval endpoints and "extrapolate" the values at those points if
%        needed.
%
%      sampleTest              - Test accuracy at arbitrary point.
%       [true]
%        false
%
%        If true, the tech should check an arbitrary point for accuracy to
%        ensure that behavior hasn't been missed, e.g., due to undersampling.
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
%   breakpoint detection, a splitting degree of 256, and a custom CHEBTECH
%   refinement function:
%      p.enableBreakpointDetection = true;
%      p.breakpointPrefs.splitLength = 256;
%      p.techPrefs.refinementFunction = @custom;
%      pref = chebpref(p);
%
%   Same thing with a slightly shorter syntax:
%      p.enableBreakpointDetection = true;
%      p.breakpointPrefs.splitLength = 256;
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
% to preferences in the upper layers (i.e., CHEBTECH, FUN, and SINGFUN) which
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

    % See above for documentation.
    properties
        maxTotalLength
        enableBreakpointDetection
        breakpointPrefs
        domain
        enableSingularityDetection
        singPrefs
        enableFunqui
        tech
        techPrefs
    end

    methods

        function p = chebpref(q)
            if ( (nargin == 1) && isa(q, 'chebpref') )
                p = q;
                return
            elseif ( nargin < 1 )
                q = struct();
            end

            % Initialize default preference values.
            p.maxTotalLength = 65536;
            p.enableBreakpointDetection = false;
                p.breakpointPrefs.splitMaxLength = 128;
                p.breakpointPrefs.splitMaxTotalLength = 6000;
            p.domain = [-1 1];
            p.enableSingularityDetection = false;
                p.singPrefs.exponentTol = 1.1*1e-11;
                p.singPrefs.maxPoleOrder = 20;
            p.enableFunqui = false;
            p.tech = 'chebtech';
            p.techPrefs = struct();
                p.techPrefs.eps = 2^(-52);
                p.techPrefs.maxLength = 65537;
                p.techPrefs.exactLength = NaN;
                p.techPrefs.extrapolate = false;
                p.techPrefs.sampleTest = true;

            % Copy fields from q, placing unknown ones in techPrefs and merging
            % incomplete substructures.
            for (field = fieldnames(q).')
                if ( isprop(p, field{1}) )
                    if ( isstruct(p.(field{1})) )
                        p.(field{1}) = chebpref.mergePrefs(p.(field{1}), ...
                            q.(field{1}));
                    else
                        p.(field{1}) = q.(field{1});
                    end
                else
                    p.techPrefs.(field{1}) = q.(field{1});
                end
            end
        end

        function out = subsref(p, ind)
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
                    if ( isprop(p, ind(1).subs) )
                        out = p.(ind(1).subs);
                    else
                        out = p.techPrefs.(ind(1).subs);
                    end

                    if ( numel(ind) > 1 )
                        out = subsref(out, ind(2:end));
                    end
                otherwise
                    error('CHEBTECH:subsref:badType', ...
                        'Invalid subscripted reference type.')
            end
        end

        function p = subsasgn(p, ind, val)
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
                    if ( isprop(p, ind(1).subs) )
                        p = builtin('subsasgn', p, ind, val);
                    else
                        p.techPrefs = builtin('subsasgn', p.techPrefs, ind, ...
                            val);
                    end
                otherwise
                    error('CHEBTECH:subsasgn:badType', ...
                        'Invalid subscripted assignment type.')
            end
        end

    end

    methods ( Static = true )

        function p = mergePrefs(p, q, map)
        %MERGEPREFS   Merge preference structures.
        %   P = MERGEPREFS(P, Q), where P and Q are MATLAB structures, "merges"
        %   Q into P by replacing the contents of fields in P with those of
        %   identically-named fields in Q.  If Q has a field whose name does
        %   not match any of those in P, it is added to P.
        %
        %   P = MERGEPREFS(P, Q, MAP) does the same but uses the structure MAP
        %   to "translate" the names of fields of Q into names of fields of P.
        %   If Q has a field FIELD and MAP has a field of the same name, then
        %   the value of P.(MAP.FIELD) will be replaced by the contents of
        %   Q.FIELD.  If P does not have a field matching the string stored in
        %   MAP.FIELD, one will be added to P.
        %
        %   P and Q may also be CHEBPREF objects.  In this case, P and Q are
        %   replaced by P.TECHPREFS and Q.TECHPREFS before proceeding.

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Developer notes:
        %  - This function is a helper function intended for use by "technology"
        %    objects (usually subclasses of SMOOTHFUN) for managing their
        %    preferences.  See CHEBTECH.PREF for an illustration.
        %  - The second syntax is useful, e.g., if Q contains abstractly-named
        %    preferences which may have a better name within the specific
        %    context of the tech object whose preferences are stored in P.
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            if ( isa(p, 'chebpref') )
                p = p.techPrefs;
            end

            if ( isa(q, 'chebpref') )
                q = q.techPrefs;
            end

            if ( nargin < 3 )
                map = struct();
            end

            for (field = fieldnames(q).')
                if ( isfield(map, field{1}) )
                    p.(map.(field{1})) = q.(field{1});
                else
                    p.(field{1}) = q.(field{1});
                end
            end
        end

    end
end

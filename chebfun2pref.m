classdef chebfun2pref < chebpref
%CHEBFUN2PREF   Class for managing preferences for Chebfun2.
%
% Constructor inputs:
%   P = chebfun2pref() creates a chebfun2pref object with the default values of the
%   preferences.  For a list of all available preferences, see above.
%
%   P = chebfun2pref(Q), where Q is a MATLAB structure uses the field/value pairs
%   in Q to set the properties of P.  If a field of Q has a name which matches
%   a property of P, the value of that property of P is set to the value
%   associated to that field in Q.  If a field of Q has a name that does not
%   correspond to a known preference, then an error is thrown.
%
%   P = chebfun2pref(Q), where Q is a chebfun2pref, sets P to be a copy of Q.
%
% See also chebfun2pref.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% TODO:  Further documentation of chebfun2pref preferences.

    methods

        function outPref = chebfun2pref(inPref)
            if ( (nargin == 1) && isa(inPref, 'chebfun2pref') )
                outPref = inPref;
                return
            elseif ( nargin < 1 )
                inPref = struct();
            end

            % Initialize default preference values.
            outPref.prefList = chebfun2pref.manageDefaultPrefs('get');

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
                    error('chebfun2pref:chebfun2pref:badPref', ...
                        'Unrecognized preference name.');
                end
            end
        end

       function display(pref)
       %DISPLAY   Display a chebfun2pref object.
       %   DISPLAY(PREF) prints out a list of the preferences stored in the
       %   chebfun2pref object PREF.

            % Compute the screen column in which pref values start.
            valueCol = 24; % length('    enableSingularityDetection:   ');

            % A subfunction to pad strings for formatting.
            function s = padString(s)
            %PADSTRING   Add whitespace to string for formatting.
                s = [s repmat(' ', 1, valueCol - length(s))];
            end

            % Print values of "known" preferences.
            prefList = pref.prefList;

            fprintf('chebfun2pref object with the following preferences:\n');
            fprintf([padString('    maxRank:') '%d\n'], ...
                prefList.maxRank');
            fprintf([padString('    maxLength:') '%d\n'], ...
                prefList.maxLength');            
            fprintf([padString('    eps:') '%d\n'], ...
                prefList.eps');            
            fprintf([padString('    exactLength:') '%d\n'], ...
                prefList.exactLength');            
            fprintf([padString('    sampleTest:') '%d\n'], ...
                prefList.sampleTest'); 
        end

    end

    methods ( Static = true )
        function pref = getFactoryDefaults(getFactory)
        %GETFACTORYDEFAULTS   Get factory default preferences.
        %   PREF = chebfun2pref.GETFACTORYDEFAULTS() returns a chebfun2pref
        %   object with the preferences set to their factory defaults,
        %   irrespective of the currently defined values of the default
        %   preferences.  This function is useful if the user wishes to
        %   solve ODEs with CHEBOP using the factory defaults when other
        %   user-set defaults are currently in force.
        %
        % See also GETDEFAULTS, SETDEFAULTS.

            fd = chebfun2pref.factoryDefaultPrefs();
            pref = chebfun2pref(fd);
        end

        function pref = getDefaults()
        %GETDEFAULTS   Get default preferences.
        %   PREF = chebfun2pref.GETDEFAULTS() returns a chebfun2pref object with
        %   the preferences set to the currently stored default values.  It is
        %   equivalent to PREF = chebfun2pref().
        %
        % See also GETFACTORYDEFAULTS, SETDEFAULTS.

            pref = chebfun2pref();
        end

        function setDefaults(varargin)
        %SETDEFAULTS   Set default preferences.
        %   chebfun2pref.SETDEFAULTS(PREF1, VAL1, PREF2, VAL2, ...) sets the
        %   default values for the preferences whose names are stored in the
        %   strings PREF1, PREF2, ..., etc. to VAL1, VAL2, ..., etc.  All
        %   subsequently constructed chebfun2pref objects will use these values
        %   as the defaults.
        %
        %   chebfun2pref.SETDEFAULTS(PREF) sets the default values to the
        %   preferences stored in the chebfun2pref object PREF.  PREF can also
        %   be a MATLAB structure, in which case it is converted to a
        %   chebfun2pref as described in the documentation for the chebfun2pref
        %   constructor first.
        %
        %   chebfun2pref.SETDEFAULTS('factory') resets the default preferences to
        %   their factory values.
        %
        % See also GETDEFAULTS, GETFACTORYDEFAULTS.

            % The reason we don't just use manageDefaults as the second
            % argument to chebpref.setDefaults and wrap it in an additional
            % anonymous function instead is to get around what seems to be a
            % bug in MATLAB.  See commit messages for more information.
            manageDefaults = @chebfun2pref.manageDefaultPrefs;
            chebpref.setDefaults(@(inPref) chebfun2pref(inPref), ...
                @(varargin) manageDefaults(varargin{:}), varargin{:});
        end
    end

    methods ( Static = true, Access = private )

        function varargout = manageDefaultPrefs(varargin)
        %MANAGEDEFAULTPREFS   Private method for handling default preferences.
        %   chebfun2pref.MANAGEDEFAULTPREFS('get') returns a structure suitable
        %   for storing in the prefList property of a chebfun2pref with all of
        %   the currently stored default preferences suitable for initializing
        %   a chebfun2pref object.
        %
        %   chebfun2pref.MANAGEDEFAULTPREFS('set-factory') restores the default
        %   preferences to their "factory" values.
        %
        %   chebfun2pref.MANAGEDEFAULTPREFS('set', PREFLIST) sets the default
        %   values to those stored in the structure PREFLIST.  PREFLIST should
        %   be a structure suitable for use as a chebfun2pref prefList.
        %
        %   chebfun2pref.MANAGEDEFAULTPREFS('set', PREF1, VAL1, PREF2, VAL2, ...)
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
        %    chebfun2pref as a MATLAB structure, and that's not something with
        %    which anyone outside of chebfun2pref should be concerned.
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            persistent defaultPrefs;

            if ( isempty(defaultPrefs) )
                defaultPrefs = chebfun2pref.factoryDefaultPrefs();
            end

            if ( strcmp(varargin{1}, 'get') )
                varargout{1} = defaultPrefs;
            elseif ( strcmp(varargin{1}, 'set-factory') )
                defaultPrefs = chebfun2pref.factoryDefaultPrefs();
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
                            error('chebfun2pref:manageDefaultPrefs:badPref', ...
                                'Unrecognized preference name.');
                        end
                        varargin(1:2) = [];
                    end
                end
            end
        end

        function factoryPrefs = factoryDefaultPrefs()
        %FACTORYDEFAULTPREFS   Get structure of factory default preferences.
        %   S = chebfun2pref.FACTORYDEFAULTPREFS() returns a structure suitable
        %   for storing in the prefList property of a chebfun2pref object that
        %   contains all of the "factory default" values of the CHEBOP
        %   preferences.

            factoryPrefs.domain = [-1 1 -1 1];
            factoryPrefs.tech = @chebtech2;
            factoryPrefs.maxRank = 513;
            factoryPrefs.maxLength = 65537;
            factoryPrefs.eps = eps;
            factoryPrefs.exactLength = 0;
            factoryPrefs.sampleTest = 1;
        end

    end

end
classdef chebpref
%CHEBPREF   Abstract class for Chebfun preferences.
%
%   CHEBPREF is an abstract class for managing preferences of the various
%   subsystems comprising Chebfun.  On its own, it provides little more than a
%   wrapper around a MATLAB struct.  Concrete preference classes inherit from
%   CHEBPREF and add to its functionality as desired.
%
% See also CHEBFUNPREF, CHEBOPPREF.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.


    properties ( Access = protected )
        % MATLAB struct to hold a list of preferences for a given subsystem.
        prefList
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Non-Static methods implemented by this class.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods ( Access = public, Static = false )

        function out = subsref(pref, ind)
        %SUBSREF   Subscripted referencing for CHEBPREF.
        %   P.PROP, where P is a CHEBPREF object, returns the value of the
        %   CHEBPREF property PROP stored in P.  If PROP is not a CHEBPREF
        %   property, an error will be thrown.
        %
        %   CHEBPREF does not support any other subscripted referencing
        %   types, including '()' and '{}'.

            switch ( ind(1).type )
                case '.'
                    if ( isfield(pref.prefList, ind(1).subs) )
                        out = pref.prefList.(ind(1).subs);
                    end

                    if ( numel(ind) > 1 )
                        out = subsref(out, ind(2:end));
                    end
                otherwise
                    error('CHEBPREF:subsref:badType', ...
                        'Invalid subscripted reference type.')
            end
        end

        function pref = subsasgn(pref, ind, val)
        %SUBSASGN   Subscripted assignment for CHEBPREF.
        %   P.PROP = VAL, where P is a CHEBPREF object, assigns the value
        %   VAL to the CHEBPREF property PROP stored in P.  If PROP is not a
        %   CHEBPREF property, an error will be thrown.
        %
        %   CHEBPREF does not support any other subscripted assignment types,
        %   including '()' and '{}'.

            switch ( ind(1).type )
                case '.'
                    if ( isfield(pref.prefList, ind(1).subs) )
                        pref.prefList = builtin('subsasgn', pref.prefList, ...
                            ind, val);
                    end
                otherwise
                    error('CHEBPREF:subsasgn:badType', ...
                        'Invalid subscripted assignment type.')
            end
        end

    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Abstract methods
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods ( Abstract = true )
        % Display information about a CHEBPREF object.
        display(pref)
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Static methods
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods ( Access = public, Static = true )

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

    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Hidden static methods.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods ( Access = public, Static = true, Hidden = true  )

        function setDefaults(makeSC, manageSCDefaults, varargin)
        %SETDEFAULTS   Set default preferences.
        %   CHEBPREF.SETDEFAULTS(MAKESC, MANAGESCDEFAULTS, ...) sets the
        %   persistently stored defaults of a subclass of CHEBPREF.  MAKESC is
        %   a handle to the subclass constructor, and MANAGESCDEFAULTS is a
        %   handle to the subclass static method responsible for managing the
        %   stored defaults.
        %
        %   Users should not call CHEBPREF.SETDEFAULTS.  Call the SETDEFAULTS
        %   method of the appropriate Chebfun subsystem whose defaults need to
        %   be altered instead.

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Developer notes:
        %  - This is a template method to prevent subclasses from needing to
        %    duplicate the logic required to set persistently-stored defaults
        %    As the help text says, this method is not user-facing.
        %  - This method is not likely to work if the subclass does not handle
        %    things in the exact same way that the existing CHEBFUNPREF and
        %    CHEBOPPREF subclasses do.
        %  - The function referred to by the MANAGESCDEFAULTS input needs to
        %    conform to a specific prototype.  See the documentation for
        %    CHEBFUNPREF.MANAGEDEFAULTPREFS for an example.
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            nargs = length(varargin);

            if ( nargs < 1)
                error('CHEBPREF:setDefaults:notEnoughArguments', ...
                    'Not enough arguments.');
            end

            if ( nargs == 1 )
                if ( isstruct(varargin{1}) )
                    varargin{1} = makeSC(varargin{1});
                end

                if ( ischar(varargin{1}) && strcmp(varargin{1}, 'factory') )
                    manageSCDefaults('set-factory');
                elseif ( isa(varargin{1}, 'chebpref') )
                    manageSCDefaults('set', varargin{1}.prefList);
                else
                    error('CHEBPREF:setDefaults:badArg', ...
                        ['When calling chebpref.setDefaults() with just ' ...
                         'one argument, that argument must be ''factory'', ' ...
                         'a subclass of CHEBPREF, or a MATLAB structure.']);
                end
            elseif ( mod(nargs, 2) == 0 )
                manageSCDefaults('set', varargin{:});
            else
                error('CHEBPREF:setDefaults:unpairedArg', ...
                    'Unpaired argument in name-value pair list.');
            end
        end

    end

end

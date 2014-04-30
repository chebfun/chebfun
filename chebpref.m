classdef chebpref
    properties ( Access = protected )
        prefList
    end

    methods

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

    methods ( Abstract = true )
        display(pref)
    end

    methods ( Static = true )

        function pref1 = mergePrefs(pref1, pref2, map)
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

    methods ( Static = true )

        function setDefaults(makeSC, manageSCDefaults, varargin)
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

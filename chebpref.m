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

end

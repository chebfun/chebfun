classdef cheboppref < chebpref
% TODO: Add CHEBFUNPREF style documentation of the options available. To be written
% when we introduce nonlinear ODEs, since those involve a number of options.
    
% TODO: The relationship between CHEBOPPREF and CHEBFUNPREF needs some serious
% consideration.

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
       %DISPLAY   Display a CHEBFUNPREF object.
       %   DISPLAY(PREF) prints out a list of the preferences stored in the
       %   CHEBFUNPREF object PREF.

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
            fd = cheboppref.factoryDefaultPrefs();
            pref = cheboppref(fd);
        end

        function pref = getDefaults()
            pref = cheboppref();
        end

        function setDefaults(varargin)
            chebpref.setDefaults(@(inPref) cheboppref(inPref), ...
                @cheboppref.manageDefaultPrefs, varargin{:});
        end
    end

    methods ( Static = true, Access = private )

        function varargout = manageDefaultPrefs(varargin)
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

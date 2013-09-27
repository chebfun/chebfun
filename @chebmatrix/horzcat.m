function C = horzcat(varargin)

for k = 1:numel(varargin)
    if ( isa(varargin{k}, 'chebmatrix') )
        varargin{k} = varargin{k}.blocks;
    end
end

C = chebmatrix( horzcat(varargin{:}) );

end
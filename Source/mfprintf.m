function mfprintf(fids,varargin)
% fprintf to multiple sources

if isempty(fids)
    return
end

for i = 1:length(fids)
    fprintf(fids(i),varargin{:});
end
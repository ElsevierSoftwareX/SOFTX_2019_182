function boundaries = extractBoundaryMulti3D_new(labels)
% Extract boundaries from max flow code (the minimum cost set). Boundary is
% the first zero voxel value along each A-scan. (this might be opposite of
% what it should be)

sz = size(labels);
if length(sz) == 3
    sz(4) = 1;
end
if length(sz) == 2
    sz(3) = 1;
    sz(4) = 1;
end

labels = labels < 1;
boundaries = zeros(sz(2:4));
for i = 1:sz(4)
    % Use the max function to get the indices of the first zero value. We
    % are relying on the max function to give the index to the ~first~ 
    % true value.
    [~,inds] = max(labels(:,:,:,i),[],1);
    
    if any(inds(:) == 1)
        npt = sum(inds(:) == 1);
        warning(sprintf('%d boundary points not found in layer %d',npt,i))
        inds(inds == 1) = nan;
    end
    
    boundaries(:,:,i) = squeeze(inds);
end

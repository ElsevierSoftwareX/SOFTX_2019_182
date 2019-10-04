function slo = octfundus(img_vol,rpe_bds)

ef_th = 20;
rpe_bds = round(rpe_bds);
ef_mask = false(size(img_vol));
for i = 1:size(img_vol,2)
    for j = 1:size(img_vol,3)
        ef_mask((rpe_bds(i,j,1)-ef_th):rpe_bds(i,j,1),i,j) = true;
    end
end
img_vol(~ef_mask) = 0;
slo = flipud(squeeze(sum(img_vol,1))');
slo = bsxfun(@rdivide,slo,mean(slo,2));
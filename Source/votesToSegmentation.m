function bd_pts = votesToSegmentation(votes_vol,seg_alg,header,f_cen)

if nargin < 3
    header = [];
end
if nargin < 4
    f_cen = [];
end

graph_file_folder = '';

if ~isempty(header)
    % The scale of the original data, to be used when calculating smoothing
    % parameters
    scale_orig = [3.87167 6.07151 129.44131];
    scale_ratio = scale_orig./...
                    (1000*[header.ScaleZ header.ScaleX header.Distance]);
    if all(abs(scale_ratio-1) < 0.3)
        % If the scale is close to the same, probably spectralis data, so don't
        % change anything (for consistency with BMOE paper work, with a fixed
        % sigma for all data)
        scale_ratio = [];
    end
else
    scale_ratio = [];
end

votes_vol = double(votes_vol);

% Convert to probabilities
if max(votes_vol(:)) > 1
    votes_vol = votes_vol/max(votes_vol(:));
end

switch seg_alg
    case 'canny'
        bd_pts = zeros(size(squeeze(votes_vol(1,:,:,:))));
        for i = 1:size(votes_vol,3)
            b = calcEdgesFromProb(squeeze(votes_vol(:,:,i,:)),[]);
            for j = 1:9
                bd_pts(:,i,j) = b{j}(:,2);
            end
        end               
    case 'optimalsurface'
        % Parameter set for graph
        % Based on fixed smoothness values used in the paper
        graph_file = [graph_file_folder 'graph_params_paper.mat'];
        
        bd_pts = optimal_surface(votes_vol,scale_ratio,graph_file);
    case 'optimalsurface_bmoe'
        % Based on fixed smoothness values used in the paper
        graph_file = [graph_file_folder 'graph_params_paper.mat'];
        
        bd_pts = optimal_surface(votes_vol,scale_ratio,graph_file);
    case 'optimalsurface_constrained'        
        % Based on maximum and minimum values from 35 subjects
        graph_file = [graph_file_folder 'graph_params_pixelwise.mat'];
        
        bd_pts = optimal_surface(votes_vol,scale_ratio,graph_file,f_cen);
end

function bd_pts = optimal_surface(votes_vol,scale_ratio,graph_file,f_cen)

if nargin < 4
    f_cen = [];
end

load(graph_file)

rs = false;
if ~isempty(scale_ratio)
    % Resize volume so it doesnt run out of memory    
    if scale_ratio(1) > 1.5
        rs = true;
        scale_ratio(1) = scale_ratio(1)/2;        
        votes_vol2 = imresize(votes_vol,[0.5*size(votes_vol,1) size(votes_vol,2)]);
        votes_vol = reshape(votes_vol2,...
            [size(votes_vol2,1) size(votes_vol2,2)...
             size(votes_vol,3) size(votes_vol,4)]);
        clear votes_vol2
    end

    % Multiplication factors due to rescaling surfaces
    scf = 1./[scale_ratio(2) scale_ratio(2) scale_ratio(3)...
           scale_ratio(3) 1 1];
    % Modify constraints for pixel size
    fields = fieldnames(params);
    for i = 1:length(fields)
        fi = params.(fields{i});
        fi = scf(i)*imresize(fi,[size(votes_vol,2) size(votes_vol,3)]);
        fi = round(fi*scale_ratio(1));
        params.(fields{i}) = fi;
    end
else
    % Should be 1024x49
    if size(votes_vol,2) ~= 1024
        % Shouldn't happen
        error('error: invalid number of A-scans, resize not yet implemented')
    end
    if size(votes_vol,3) ~= 49
        % Add or remove slices
        sdiff = size(votes_vol,3)-49;
        sdiff1 = abs(ceil(sdiff/2));
        sdiff2 = abs(floor(sdiff/2));
        fields = fieldnames(params);
        for i = 1:length(fields)
            fi = params.(fields{i});
            if sdiff > 0
                % before first image
                fi = cat(2,repmat(fi(:,1,:),[1 sdiff1 1]),fi);
                % after last image
                fi = cat(2,fi,repmat(fi(:,end,:),[1 sdiff2 1]));
            else
                fi = fi(:,(1+sdiff1):(end-sdiff2),:);
            end
            params.(fields{i}) = fi;
        end
    end
end

if ~isempty(f_cen)
    % Center constraints at the fovea
    vs = fliplr(size(params.DxL(:,:,1)));
    
    % Align votes to center of graph
    g_cen = ceil(vs/2);
    
    cen_diff = g_cen - f_cen;
    
    fields = fieldnames(params);
    for i = 1:length(fields)
        fi = params.(fields{i});
        if cen_diff(1) > 0
            % Add to end, remove from beginning
            fi = cat(2,fi,repmat(fi(:,end,:),[1 cen_diff(1) 1]));
            fi(:,1:cen_diff(1),:) = [];
        elseif cen_diff(1) < 0
            % Add to beginning, remove from end
            fi = cat(2,repmat(fi(:,1,:),[1 -cen_diff(1) 1]),fi);
            fi(:,(end+cen_diff(1)+1):end,:) = [];
        end

        if cen_diff(2) > 0
            % Add to end, remove from beginning
            fi = cat(1,fi,repmat(fi(end,:,:),[cen_diff(2) 1 1]));
            fi(1:cen_diff(2),:,:) = [];
        elseif cen_diff(2) < 0
            % Add to beginning, remove from end
            fi = cat(1,repmat(fi(1,:,:),[-cen_diff(2) 1 1]),fi);
            fi((end+cen_diff(2)+1):end,:,:) = [];
        end
        
        params.(fields{i}) = fi;
    end
end

% Allow 0 thickness layers
dl = params.dl;
dl(dl == 1) = 0;
params.dl = dl;

if ~isempty(scale_ratio) && rs
    bd_pts = 2*calcBoundariesFromOptimalSurface(votes_vol,params);
else
    bd_pts = calcBoundariesFromOptimalSurface(votes_vol,params);
end
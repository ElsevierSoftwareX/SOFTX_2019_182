function pred_pts = calcEdgesFromProb(prob_img,img)

pred_pts = cell(1,size(prob_img,3)-1);

% Parameters
% -Sigma for gaussian smoothing of probabilities
sigma = 1;
% -Hysteresis lower and upper threshold (percentage of min and max
%   intensities)
lower_ratio = 0.1;
upper_ratio = 0.5;

prob_img = imfilter(prob_img,fspecial('gaussian',2*round(3*sigma)+1,sigma),'same','replicate');

layerOrder = [1 2 9 8 7 6 3 5 4];
for i = layerOrder    
    layer_prob = prob_img(:,:,i);        

    % Non-maximal suppression
    nms_img = nms(layer_prob);        

    upper = upper_ratio*max(nms_img(:));
    lower = lower_ratio*max(nms_img(:));

    % Hysteresis
    hyst_img = hysthresh(nms_img,lower,upper);

    % Convert hysteresis threshold to have a single point per ascan
    pts = hystToPoints(hyst_img);

    % Interpolate nan values
    pts = interpVals(pts,'spline');


    % Enforce topology of each a-scan
    if i < size(prob_img,3)
        % From layers below
        f = find(~cellfun(@isempty,pred_pts(i+1:end)),1,'first');
        if ~isempty(f)
            pts_post = pred_pts{i+f};
            pts(pts(:,2) >= pts_post(:,2),2) = pts_post(pts(:,2) >= pts_post(:,2),2) - 1;
        end
        % From layers above
        f = find(~cellfun(@isempty,pred_pts(1:i-1)),1,'last');
        if ~isempty(f)
            pts_prev = pred_pts{f};
            pts(pts(:,2) <= pts_prev(:,2),2) = pts_prev(pts(:,2) <= pts_prev(:,2),2) + 1;
        end
    end

    pred_pts{i} = pts;
end


function pts = interpVals(pts,method)

% Interpolate nan values
if any(isnan(pts(:,2)))
nanpts = find(isnan(pts(:,2)));

% Find nan points at the beginning
startflag = true; ni = 1;
while startflag
    if (length(nanpts)>=ni) && (nanpts(ni) == ni)
        ni = ni + 1;
    else 
        startflag = false;
    end                
end
% Send equal to first non-nan point
pts(1:(ni-1),2) = pts(ni,2);
nanpts(1:(ni-1)) = [];

% Find nan points at the end
endflag = true; ni = 1;
while endflag
    if (length(nanpts)>=ni) && (nanpts(end-ni+1) == size(pts,1)-ni+1);                    
        ni = ni + 1;
    else 
        endflag = false;
    end                
end
% Send equal to last non-nan point
pts((end-ni+2):end,2) = pts(end-ni+1,2);
nanpts((end-ni+2):end) = [];

if ~isempty(nanpts)
    % Interpolate with one knot every 30 points
    subsamp = 30;
    np = isnan(pts(:,2));

    % Look for start and end of nan sequences, make sure we use the
    % last non-nan point before and first non-nan point after for our
    % interpolation 
    dp = diff(np);
    sp = find(dp==1); ep = find(dp==-1)+1;

    % Subsample points
    fp = find(~np);
    fp = unique([fp(1:subsamp:end); sp; ep]); % unique sorts

    pts(nanpts,2) = round(interp1(pts(fp,1),pts(fp,2),pts(np,1),method));
end            

end

function pts = hystToPoints(hyst_img)

% Algorithm:    
%   - Traverse each ascan
%   - If we get to an ascan with multiple points
%       - Choose point closest to previous point
%   - If we get to an ascan with only one point and a big jump
%       - Check if the previous ascan had multiple points and
%         retroactivly change it if we made a wrong decision before
pts = [1:size(hyst_img,2); zeros(1,size(hyst_img,2))]';
max_diff = 3;
j = 1;
bt_ind = 0;
while j <= size(hyst_img,2)
    inds = find(hyst_img(:,j)>0);
    if ~isempty(inds)
        if size(inds,1) > 1
            % Multiple points to choose from
            if (j == 1)
                % At the beginning, choose the point with the highest
                %   probability
                [~, indmax] = max(hyst_img(inds,j));
                pts(j,2) = inds(indmax);                
            elseif isnan(pts(j-1,2))
                % No previous point, pick the closest previous point
                lp = find(~isnan(pts(1:(j-1),2)),1,'last');
                if isempty(lp)
                    [~, indmax] = max(hyst_img(inds,j));
                    pts(j,2) = inds(indmax);
                else
                    prev_pt = pts(lp,2);
                    diff_pt = abs(inds - prev_pt);
                    [minval, indmin] = min(diff_pt);
                    if minval < max_diff
                        pts(j,2) = inds(indmin);
                    else
                        % Too far away, choose highest probability
                        [~, indmax] = max(hyst_img(inds,j));
                        pts(j,2) = inds(indmax);
                    end
                end
            else
                % Look at the previous point and choose the current point
                %   closest to the previous
                prev_pt = pts(j-1,2);
                diff_pt = abs(inds - prev_pt);

                [minval, indmin] = min(diff_pt);
                if minval < max_diff || any(j == bt_ind)
                    % We're close enough to the previous point, add it
                    pts(j,2) = inds(indmin);
                else
                    % Backtrack! see below - note that it is now a
                    % combinatorial problem with multiple points at
                    % the current index and possibly multiple
                    % points at the previous index, to temporarily
                    % combat this, choose the closest point to the
                    % previous point, the hope being that if this
                    % is a bad point, we will back track in the
                    % future to correct it
                    bt_ind(end+1) = j;
                    pts(j,2) = inds(indmin);
                    cur_ind = inds(indmin);
                    while j > 1 %&& abs(cur_ind - pts(j-1,2)) >= max_diff
                        inds_prev = find(hyst_img(:,j-1)>0);
                        if (j == 1) || ~isempty(inds_prev)
                            diff_pt = abs(inds_prev - cur_ind);
                            [minval, indmin] = min(diff_pt);
                            if minval >= max_diff
                                % We are still too far away, lets just go
                                %   back to where we were
                                j = j - 1;
                                break
                            else
                                if (inds_prev(indmin) ~= pts(j-1,2))
                                    % Backtrack!
                                    pts(j-1,2) = inds_prev(indmin);
                                    cur_ind = inds_prev(indmin);
                                    j = j - 1;
                                else
                                    if abs(cur_ind - pts(j-1,2)) >= max_diff
                                        % We back tracked and it didn't help, lets
                                        % go forward again (redo what we already
                                        % did)
                                        j = j - 1;
                                    else
                                        % Good job! Lets go back to where we were
                                        %   before
                                        j = bt_ind(end);
                                    end
                                    break
                                end
                            end
                        else
                            % We hit the wall, lets hope we're good and go
                            %   back to where we were!
                            j = bt_ind(end);
                            break
                        end
                    end
                end
            end
        else
            % Only one point to choose from
            if (j == 1) || isnan(pts(j-1,2))
                % We are either at the beginning or there is no previous
                %   point so we have to choose this point (although we
                %   might back track later)
                pts(j,2) = inds;
            else 
                % Compare to previous point
                if abs(inds - pts(j-1,2)) < max_diff || any(j == bt_ind)
                    % We are close enough to the previous point so accept
                    %   it!
                    pts(j,2) = inds;
                else
                    bt_ind(end+1) = j;
                    
                    % Lets first try to push forward for a few points and
                    % see if we can get back on track. The missing points
                    % in between can be interpolated.
                    jf = 1;
                    pth = false;
                    while (j+jf) < size(hyst_img,2) && (jf < 10)
                        inds_next = find(hyst_img(:,j+jf)>0);
                        diff_pt = abs(inds_next - pts(j-1,2));
                        [minval, indmin] = min(diff_pt);
                        if minval < max_diff
                            % Success! 
                            bt_ind(end) = [];
                            pts(j:(j+jf-1),2) = nan;
                            pth = true; % Pushed through!
                            j = j + jf -1;
                            break;
                        end 
                        jf = jf + 1;
                    end                        
                    
                    % Set this point as truth, but backtrack, maybe
                    %   we made a mistake previously (think fork in
                    %   the road with one direction leading to a
                    %   dead end)                    
                    while j > 1 && ~pth %&& abs(cur_ind - pts(j-1,2)) >= max_diff
                        pts(j,2) = inds;
                        cur_ind = inds;
                        inds_prev = find(hyst_img(:,j-1)>0);
                        if (j == 1) || ~isempty(inds_prev)
                            diff_pt = abs(inds_prev - cur_ind);
                            [minval, indmin] = min(diff_pt);
                            if minval >= max_diff
                                % We are still too far away, lets just go
                                %   back to where we were
                                j = j - 1;
                                break
                            else
                                if (inds_prev(indmin) ~= pts(j-1,2))
                                    % Backtrack!
                                    pts(j-1,2) = inds_prev(indmin);
                                    cur_ind = inds_prev(indmin);
                                    j = j - 1;
                                else
                                    if abs(cur_ind - pts(j-1,2)) >= max_diff
                                        % We back tracked and it didn't help, lets
                                        % go forward again (redo what we already
                                        % did)
                                        j = j - 1;
                                    else
                                        % Good job! Lets go back to where we were
                                        %   before
                                        j = bt_ind(end);
                                    end
                                    break
                                end
                            end
                        else
                            % We hit the wall, lets hope we're good and go
                            %   back to where we were!
                            j = bt_ind(end);
                            break
                        end
                    end
                end
            end
        end
    else
        % Did not find a point, interpolate it later
        pts(j,2) = nan;
    end
    j = j + 1;
end


function nms_img = nms(img)
% Non-maximal supression along each a-scan (vertical direction only)
%   Not memory efficient but fast, works in 3D

% Allocate memory
nms_img = zeros(size(img));

% Create volumes for one pixel above and below 
padvals = repmat(img(1,:),[1 1]);        
above_p = cat(1,padvals,img(1:end-1,:));
padvals = repmat(img(end,:,:),[1 1 1]);
below_p = cat(1,img(2:end,:,:),padvals);  

% Compare to pixels above and below
nms_img(img >= above_p & img >= below_p) = 1;
nms_img = nms_img.*img;

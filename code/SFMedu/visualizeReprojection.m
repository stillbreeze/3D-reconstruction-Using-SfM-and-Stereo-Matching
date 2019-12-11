function visualizeReprojection(graph, frames)

nfigs = size(graph.ObsIdx,1);
w = ceil(sqrt(nfigs));
h = floor(sqrt(nfigs));
subplot(w,h,1);
for c=1:size(graph.ObsIdx,1)
    subplot(w,h,c);
    im = imresize(imread(frames.images{c}),frames.imsize(1:2));
    imshow(im);
    hold on
    X = f2K(graph.f) * transformPtsByRt(graph.Str,graph.Mot(:,:,c));
    xy = X(1:2,:) ./ X([3 3],:);
    selector = find(graph.ObsIdx(c,:)~=0);
    unselector = graph.ObsIdx(c,:)==0;
    tx_points = xy(:, selector);
    unmatching_points = xy(:, unselector);
    
    matching_points = graph.ObsVal(:,graph.ObsIdx(c,selector));
    fprintf('matching points: %d\n', size(matching_points));
    
    centered_tx_points(1, :) = size(im, 2)/2 - tx_points(1, :);
    centered_tx_points(2, :) = size(im, 1)/2 - tx_points(2, :);
    centered_matching_points(1, :) = size(im, 2)/2 - matching_points(1, :);
    centered_matching_points(2, :) = size(im, 1)/2 - matching_points(2, :);
    centered_unmatching_points(1, :) = size(im, 2)/2 - unmatching_points(1, :);
    centered_unmatching_points(2, :) = size(im, 1)/2 - unmatching_points(2, :);
    
    plot(centered_tx_points(1, :),centered_tx_points(2, :), 'g+');    
    plot(centered_matching_points(1, :), centered_matching_points(2, :), 'rx');
    linesx = [centered_tx_points(1, :); centered_matching_points(1, :)];
    linesy = [centered_tx_points(2, :); centered_matching_points(2, :)];
    plot(linesx, linesy, 'g');
    plot(centered_unmatching_points(1, :), centered_unmatching_points(2, :), 'yo');
    
    clear centered_tx_points;
    clear centered_matching_points;
    clear centered_unmatching_points;
end


end
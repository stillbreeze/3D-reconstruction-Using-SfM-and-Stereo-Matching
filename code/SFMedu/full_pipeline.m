function full_pipeline(data_seq_idx, adjust_focal_length, do_sequential_match, maxSize, visualize, do_dense)
%% run full pipeline with specified parameters

%% data

if data_seq_idx == 0
    frames.images{1}='images/B21.jpg';
    frames.images{2}='images/B22.jpg';
    frames.images{3}='images/B23.jpg';
    frames.images{4}='images/B24.jpg';
    frames.images{5}='images/B25.jpg';
    try
        frames.focal_length = extractFocalFromEXIF(frames.images{1});
    catch
    end
    if ~isfield(frames,'focal_length') || isempty(frames.focal_length)
        fprintf('Warning: cannot find the focal length from the EXIF\n');
        frames.focal_length = 719.5459; % for testing with the B??.jpg sequences
    end
elseif  data_seq_idx == 1
    frames.images{1}='images/couch21.jpg';
    frames.images{2}='images/couch22.jpg';
    frames.images{3}='images/couch23.jpg';
    frames.images{4}='images/couch24.jpg';
    frames.images{5}='images/couch25.jpg';
    frames.focal_length = 968.3333333333334 ;%1106.67;
% % 
elseif data_seq_idx == 2
    frames.images{1}='images/mug1.jpg';
    frames.images{2}='images/mug2.jpg';
    frames.images{3}='images/mug3.jpg';
    frames.images{4}='images/mug4.jpg';
    frames.images{5}='images/mug5.jpg';
    frames.focal_length = 968.3333333333334 ;%1106.67;
elseif data_seq_idx == 3
    frames.images{1}='images/self1_1.jpg';
    frames.images{2}='images/self1_2.jpg';
    frames.images{3}='images/self1_3.jpg';
    frames.images{4}='images/self1_4.jpg';
    frames.images{5}='images/self1_5.jpg';
    frames.focal_length = 3300 ;
elseif data_seq_idx == 4
    frames.images{1}='images/self2_1.jpg';
    frames.images{2}='images/self2_2.jpg';
    frames.images{3}='images/self2_3.jpg';
    frames.images{4}='images/self2_4.jpg';
    frames.images{5}='images/self2_5.jpg';
    frames.focal_length = 3300 ;
elseif data_seq_idx == 5
    frames.images{1}='images/self3_1.jpg';
    frames.images{2}='images/self3_2.jpg';
    frames.images{3}='images/self3_3.jpg';
    frames.images{4}='images/self3_4.jpg';
    frames.images{5}='images/self3_5.jpg';
    frames.focal_length = 3300 ;
elseif data_seq_idx == 6
    frames.images{1}='images/self4_1.jpg';
    frames.images{2}='images/self4_2.jpg';
    frames.images{3}='images/self4_3.jpg';
    frames.images{4}='images/self4_4.jpg';
    frames.images{5}='images/self4_5.jpg';
    frames.focal_length = 3300 ;
else
   error('unknown dataseq idx\n'); 
end

%% data 
frames.length = length(frames.images);



frames.imsize = size(imread(frames.images{1}));
if max(frames.imsize)>maxSize
    scale = maxSize/max(frames.imsize);
    frames.focal_length = frames.focal_length * scale;
    frames.imsize = size(imresize(imread(frames.images{1}),scale));
end


frames.K = f2K(frames.focal_length);
disp('intrinsics:');
disp(frames.K);

all_pairs = [];
adj_mat = zeros(frames.length);


%% do MST match... all pairs matching + weighted Nkeypoint minimum spanning tree
if ~do_sequential_match
    fprintf('\n\n\tDoing Minimum Spanning Tree matching \n\n');
    for frame=1:frames.length-1
        for f2=frame+1:frames.length
            s = RandStream('mcg16807','Seed',10); 
            RandStream.setGlobalStream(s);

             % higher the score, the more matches returned (easier to
             % match with given score)
            [pair, score] = match2viewSIFTReturnThresh(frames, frame, f2);
            iscore = 1 - score;
            pair.nmatches = size(pair.matches, 2);
            if ~isstruct(all_pairs)
                all_pairs = repmat(pair, frames.length + 1, frames.length + 1);
            end
            all_pairs(frame,f2) = pair;
            all_pairs(f2,frame) = pair;
            adj_mat(frame,f2) = pair.nmatches * iscore;
            adj_mat(f2,frame) = pair.nmatches * iscore;

        end
    end
    
    maxmat = max(max(adj_mat));
    adj_mat_flip = maxmat + 1 - adj_mat;
    adj_mat_flip(adj_mat_flip == maxmat + 1) = 0;
    [mst, trav] = graphminspantree(sparse(adj_mat_flip));
    path = graphpred2path(trav);
    full_path = path{length(path)};
    ordered_inds = zeros(length(full_path)-1, 2);
    for ii=1:length(full_path) - 1
       ordered_inds(ii, :) = [full_path(ii), full_path(ii+1)];
    end
    
    count = 0;
    for oidx = 1:length(ordered_inds)
        oi = ordered_inds(oidx, :);
        count = count + 1;
        % need to set this random seed to produce exact same result
        s = RandStream('mcg16807','Seed',10); 
        RandStream.setGlobalStream(s);

        fprintf('adding edge (%d, %d)', oi(1), oi(2));
        pair = all_pairs(oi(1), oi(2));

%         if visualize, showMatches(pair,frames); title('raw feature matching'); end

        if true % choose between different ways of getting E
            % Estimate Fundamental matrix
            pair = estimateF(pair);    
            % Convert Fundamental Matrix to Essential Matrix
            pair.E = frames.K' * pair.F * frames.K; % MVG Page 257 Equation 9.12
        else
            % Estimate Essential Matrix directly using 5-point algorithm
            pair = estimateE(pair,frames); 
        end


%         if visualize, showMatches(pair,frames); title('inliers'); end

        % Get Poses from Essential Matrix
        pair.Rt = RtFromE(pair,frames);

        % Convert the pair into the BA format
        Graph{count} = pair2graph(pair,frames);

        % re-triangulation
        Graph{count} = triangulate(Graph{count},frames);
%         if visualize, visualizeGraph(Graph{count},frames); title('triangulation'); end

        % outlier rejection
        % Graph{frame} = removeOutlierPts(Graph{frame});

        % bundle adjustment
        if adjust_focal_length
            Graph{count} = bundleAdjustmentFull(Graph{count});
        else
            Graph{frame} = bundleAdjustment(Graph{frame}, adjust_focal_length);
        end
%         if visualize, visualizeGraph(Graph{count},frames); title('after two-view bundle adjustment'); end
    end
end


        



%% SIFT matching and Fundamental Matrix Estimation
if do_sequential_match
    fprintf('\n\n\tDoing Sequential matching \n\n');
    for frame=1:frames.length-1    
        % need to set this random seed to produce exact same result
        s = RandStream('mcg16807','Seed',10); 
        RandStream.setGlobalStream(s);

        % keypoint matching
        %pair = match2viewSURF(frames, frame, frame+1);
        pair = match2viewSIFT(frames, frame, frame+1);

%         if visualize, showMatches(pair,frames); title('raw feature matching'); end

        if true % choose between different ways of getting E
            % Estimate Fundamental matrix
            pair = estimateF(pair);    
            % Convert Fundamental Matrix to Essential Matrix
            pair.E = frames.K' * pair.F * frames.K; % MVG Page 257 Equation 9.12
        else
            % Estimate Essential Matrix directly using 5-point algorithm
            pair = estimateE(pair,frames); 
        end


%         if visualize, showMatches(pair,frames); title('inliers'); end

        % Get Poses from Essential Matrix
        pair.Rt = RtFromE(pair,frames);

        % Convert the pair into the BA format
        Graph{frame} = pair2graph(pair,frames);

        % re-triangulation
        Graph{frame} = triangulate(Graph{frame},frames);
%         if visualize, visualizeGraph(Graph{frame},frames); title('triangulation'); end

        % outlier rejection
        % Graph{frame} = removeOutlierPts(Graph{frame});

        % bundle adjustment
        if adjust_focal_length
            Graph{frame} = bundleAdjustmentFull(Graph{frame});
        else
            Graph{frame} = bundleAdjustment(Graph{frame}, adjust_focal_length);
        end
%         if visualize, visualizeGraph(Graph{frame},frames); title('after two-view bundle adjustment'); end
    end
end

drawnow
fprintf('visualize: %d\n', visualize);

%% merge the graphs
%close all
fprintf('\n\nmerging graphs....\n');

sizes = zeros(1, length(Graph));
for i=1:length(Graph)
   sizes(i) = size(Graph{i}.matches, 2);
end

% original used for naive KP-match sorted order, but fails when the pair
% isn't computed (no edge between two frames)
% [~, sort_inds] = sort(sizes, 'descend');

% sorted order, which may be sequential, or corresponds to KP MST
sort_inds = 1:length(Graph);
mergedGraph = Graph{sort_inds(1)};

for idx=2:frames.length-1
    frame = sort_inds(idx);
    % merge graph
    mergedGraph = merge2graphs(mergedGraph,Graph{frame});
    
    % re-triangulation
    mergedGraph = triangulate(mergedGraph,frames);
%     if visualize, visualizeGraph(mergedGraph,frames); title('triangulation'); end
    
    % outlier rejection
    % mergedGraph = removeOutlierPts(mergedGraph,10);
    
    % bundle adjustment

    if adjust_focal_length
        mergedGraph = bundleAdjustmentFull(mergedGraph);
    else
        mergedGraph = bundleAdjustment(mergedGraph, adjust_focal_length);
    end
        
    
    % outlier rejection
    mergedGraph = removeOutlierPts(mergedGraph, 10);
    
    % bundle adjustment

    if adjust_focal_length
        mergedGraph = bundleAdjustmentFull(mergedGraph);
    else
        mergedGraph = bundleAdjustment(mergedGraph, adjust_focal_length);         
    end
%     if visualize, visualizeGraph(mergedGraph,frames); title('after bundle adjustment'); end
end

printReprojectionError(mergedGraph); % [for homework]

if visualize
    f = figure();
    visualizeReprojection(mergedGraph,frames);
    drawnow;
    features_file_name = char(strcat(string(data_seq_idx),'_feat.png'));
    saveas(f,features_file_name);
end% [for homework] 


points2ply('sparse.ply',mergedGraph.Str);

if frames.focal_length ~= mergedGraph.f
    disp('Focal length is adjusted by bundle adjustment');
    frames.focal_length = mergedGraph.f;
%     frames.K = f2K(frames.focal_length);
    frames.K = graph2K(mergedGraph);
    disp(frames.K);
end


%% dense matching
if do_dense
    fprintf('dense matching ...\n');
    for frame=1:frames.length-1
        Graph{frame} = denseMatch(Graph{frame}, frames, frame, frame+1);
    end


    %% dense reconstruction
    fprintf('triangulating dense points ...\n');
    for frame=1:frames.length-1
        clear X;
        P{1} = frames.K * mergedGraph.Mot(:,:,frame);
        P{2} = frames.K * mergedGraph.Mot(:,:,frame+1);
        %par
        for j=1:size(Graph{frame}.denseMatch,2)
            X(:,j) = vgg_X_from_xP_nonlin(reshape(Graph{frame}.denseMatch(1:4,j),2,2),P,repmat([frames.imsize(2);frames.imsize(1)],1,2));
        end
        X = X(1:3,:) ./ X([4 4 4],:);
        x1= P{1} * [X; ones(1,size(X,2))];
        x2= P{2} * [X; ones(1,size(X,2))];
        x1 = x1(1:2,:) ./ x1([3 3],:);
        x2 = x2(1:2,:) ./ x2([3 3],:);
        Graph{frame}.denseX = X;
        Graph{frame}.denseRepError = sum(([x1; x2] - Graph{frame}.denseMatch(1:4,:)).^2,1);

        Rt1 = mergedGraph.Mot(:, :, frame);
        Rt2 = mergedGraph.Mot(:, :, frame+1);
        C1 = - Rt1(1:3, 1:3)' * Rt1(:, 4);
        C2 = - Rt2(1:3, 1:3)' * Rt2(:, 4);
        view_dirs_1 = bsxfun(@minus, X, C1);
        view_dirs_2 = bsxfun(@minus, X, C2);
        view_dirs_1 = bsxfun(@times, view_dirs_1, 1 ./ sqrt(sum(view_dirs_1 .* view_dirs_1)));
        view_dirs_2 = bsxfun(@times, view_dirs_2, 1 ./ sqrt(sum(view_dirs_2 .* view_dirs_2)));
        Graph{frame}.cos_angles = sum(view_dirs_1 .* view_dirs_2);

        c_dir1 = Rt1(3, 1:3)';
        c_dir2 = Rt2(3, 1:3)';
        Graph{frame}.visible = (sum(bsxfun(@times, view_dirs_1, c_dir1)) > 0) & (sum(bsxfun(@times, view_dirs_2, c_dir2)) > 0);
    end

    % visualize the dense point cloud
    if visualize
        figure
        for frame=1:frames.length-1
            hold on
            goodPoint =  Graph{frame}.denseRepError < 0.05;
            plot3(Graph{frame}.denseX(1,goodPoint),Graph{frame}.denseX(2,goodPoint),Graph{frame}.denseX(3,goodPoint),'.b','Markersize',1);
        end
        hold on
        plot3(mergedGraph.Str(1,:),mergedGraph.Str(2,:),mergedGraph.Str(3,:),'.r')
        axis equal
        title('dense cloud')
        for i=1:frames.length
            drawCamera(mergedGraph.Mot(:,:,i), frames.imsize(2), frames.imsize(1), frames.K(1,1), 0.001,i*2-1);
        end
        axis tight
    end

    % output as ply file to open in Meshlab (Open Software available at http://meshlab.sourceforge.net )
    plyPoint = [];
    plyColor = [];
    for frame=1:frames.length-1
        goodPoint =  (Graph{frame}.denseRepError < 0.05) & (Graph{frame}.cos_angles < cos(5 / 180 * pi)) & Graph{frame}.visible;
        X = Graph{frame}.denseX(:,goodPoint);
        % get the color of the point
        P{1} = frames.K * mergedGraph.Mot(:,:,frame);
        x1= P{1} * [X; ones(1,size(X,2))];
        x1 = round(x1(1:2,:) ./ x1([3 3],:));
        x1(1,:) = frames.imsize(2)/2 - x1(1,:);
        x1(2,:) = frames.imsize(1)/2 - x1(2,:);
        indlin = sub2ind(frames.imsize(1:2),x1(2,:),x1(1,:));
        im = imresize(imread(frames.images{frame}),frames.imsize(1:2));
        imR = im(:,:,1);
        imG = im(:,:,2);
        imB = im(:,:,3);
        colorR = imR(indlin);
        colorG = imG(indlin);
        colorB = imB(indlin);
        plyPoint = [plyPoint X];
        plyColor = [plyColor [colorR; colorG; colorB]];
    end
    dense_file_name = strcat(string(data_seq_idx),'_dense.ply');
    points2ply(dense_file_name,plyPoint,plyColor);
    fprintf('SFMedu is finished.\n Open the result dense.ply in Meshlab (Open Software available at http://meshlab.sourceforge.net ).\n Enjoy!\n');
else
    fprintf('Not running dense reconstruction. Done\n');
end

end


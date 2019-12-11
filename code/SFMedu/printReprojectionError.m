function printReprojectionError(graph)

Error = 0;
for c=1:size(graph.ObsIdx,1)
    X = f2K(graph.f) * transformPtsByRt(graph.Str,graph.Mot(:,:,c));
    xy = X(1:2,:) ./ X([3 3],:);
    selector = find(graph.ObsIdx(c,:)~=0);
    diff = xy(:,selector) - graph.ObsVal(:,graph.ObsIdx(c,selector));
    Error = Error + norm(diff);
end
fprintf('total reprojection error: %.3f\n', E);
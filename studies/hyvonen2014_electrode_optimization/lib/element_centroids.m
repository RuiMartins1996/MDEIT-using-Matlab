function centroids = element_centroids(ctx)
%ELEMENT_CENTROIDS  [n_elem x 2] element centroids of the fixed mesh.
centroids = (ctx.nodes(ctx.elems(:,1),:) + ctx.nodes(ctx.elems(:,2),:) + ...
             ctx.nodes(ctx.elems(:,3),:)) / 3;
end

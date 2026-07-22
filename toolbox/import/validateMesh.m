function results = validateMesh(V, F)

numFaces = size(F,1);
numVertices = size(V,1);

% Basic index sanity before building anything
assert(all(F(:) >= 1) && all(F(:) <= numVertices), ...
    'Face indices out of range [1, numVertices].');

% Degenerate triangles (repeated vertex in a face -> zero area)
degenerateMask = (F(:,1) == F(:,2)) | (F(:,2) == F(:,3)) | (F(:,3) == F(:,1));
hasDegenerateFaces = any(degenerateMask);

% Directed edges per triangle: (v1->v2), (v2->v3), (v3->v1)
edgesDir = [F(:,1), F(:,2); F(:,2), F(:,3); F(:,3), F(:,1)];
adjDir = sparse(edgesDir(:,1), edgesDir(:,2), 1, numVertices, numVertices);
adjUnDir = adjDir + adjDir';

% Edge-manifold / watertight: every edge shared by exactly 2 faces
[rows, cols, vals] = find(adjUnDir);
uniqueEdgeIdx = rows < cols;
edgeCounts = vals(uniqueEdgeIdx);
edgeManifold = ~any(edgeCounts > 2);
watertight = ~any(edgeCounts ~= 2);

% Orientability: adjacent triangles must traverse shared edges in
% opposite directions, i.e. adjDir should equal its transpose's "mirror"
% with no edge duplicated in the same direction twice
D = adjDir - adjDir';
orientable = ~any(D(:)) && ~any(adjDir(:) > 1);

% Unreferenced vertices
referencedVerts = unique(F(:));
hasUnreferencedVertices = (numel(referencedVerts) ~= numVertices);

% True vertex-manifoldness: the triangle fan around each vertex must
% form a single connected loop/chain, not multiple disjoint fans meeting
% only at a point ("bowtie")
vertexManifold = checkVertexManifold(F, numVertices);

results = struct( ...
    'edgeManifold', edgeManifold, ...
    'watertight', watertight, ...
    'orientable', orientable, ...
    'vertexManifold', vertexManifold, ...
    'hasDegenerateFaces', hasDegenerateFaces, ...
    'hasUnreferencedVertices', hasUnreferencedVertices);
end


function ok = checkVertexManifold(F, numVertices)
% For each vertex, collect incident faces and check that they form a
% single connected fan (faces sharing an edge at that vertex are linked).
% A bowtie vertex has 2+ disconnected face-groups meeting only at the
% vertex itself.

    faceIdx = repmat((1:size(F,1))', 3, 1);
    vertIdx = F(:);
    [sortedV, order] = sort(vertIdx);
    sortedF = faceIdx(order);

    ok = true;
    n = numel(sortedV);
    i = 1;
    while i <= n
        j = i;
        while j <= n && sortedV(j) == sortedV(i)
            j = j + 1;
        end
        incidentFaces = sortedF(i:j-1);
        v = sortedV(i);
        if numel(incidentFaces) > 1
            if ~isFanConnected(F, incidentFaces, v)
                ok = false;
                return;
            end
        end
        i = j;
    end
end


function connected = isFanConnected(F, incidentFaces, v)
% Two incident faces are "linked" at v if they share an edge that
% contains v (i.e. share the OTHER vertex opposite v in their local
% triangle, one of the two edges at v). Build adjacency among the
% incident faces themselves and check single-component connectivity.

    m = numel(incidentFaces);
    if m <= 1
        connected = true;
        return;
    end

    % For each incident face, get the pair of "other" vertices (the
    % triangle minus v) — these are the two edge-endpoints touching v.
    others = zeros(m, 2);
    for k = 1:m
        tri = F(incidentFaces(k), :);
        tri = tri(tri ~= v);
        others(k, :) = tri;
    end

    % Build adjacency: two faces are linked if they share one of these
    % "other" vertices
    adj = false(m, m);
    for a = 1:m
        for b = a+1:m
            if any(others(a,1) == others(b,:)) || any(others(a,2) == others(b,:))
                adj(a,b) = true;
                adj(b,a) = true;
            end
        end
    end

    visited = false(m,1);
    stack = 1;
    visited(1) = true;
    while ~isempty(stack)
        cur = stack(end);
        stack(end) = [];
        neighbors = find(adj(cur,:) & ~visited');
        visited(neighbors) = true;
        stack = [stack, neighbors];
    end
    connected = all(visited);
end

function neighbors = get_neighbors(graph, rootid)
%GET_NEIGHBORS Outputs neighbors of ROOTID within GRAPH.
neighbors = {};
i = 1;
for nodeid=1:size(graph, 2)
    if (rootid ~= nodeid) && (graph(rootid, nodeid) == 1)
        neighbors{i} = nodeid;
        i = i + 1;
    end
end    
end

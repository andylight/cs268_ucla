function [comps, varargout] = connected_components(graph)
%CONNECTED_COMPONENTS Computes the connected components of an input
%graph, where GRAPH is an adjacency matrix.
%Output
%comps = connected_components(graph)
%   cell COMPS: {ids_1, ids_2, ...}
%       Where each ids_i = {id_1, id_2, ...}
%[comps, subgraphs] = connected_components(graph)
nb_nodes = size(graph, 1);
history = containers.Map('KeyType', 'int32', 'ValueType', 'logical');
unvisited = java.util.Stack();

comps = {};

for nodeid=1:nb_nodes
    unvisited.push(nodeid);
end

while ~unvisited.isEmpty()
    rootid = unvisited.pop();
    if history.isKey(rootid)
        continue;
    end
    current = java.util.Stack();
    current.push(rootid);
    cur_comp = {};
    while ~current.isEmpty()
        nodeid = current.pop();
        if history.isKey(nodeid)
            continue;
        end
        history(nodeid) = true;
        cur_comp = [cur_comp, nodeid];
        neighbors = get_neighbors(graph, nodeid);
        for i=1:length(neighbors)
            current.push(neighbors{i});
        end
    end
    comps{length(comps)+1} = cur_comp;
end

if nargout == 2
    varargout{1} = subgraphs;
end

end

function subgraphs = compute_subgraphs(graph, comps)
%COMPUTE_SUBGRAPH Computes the N graphs for each of the N connected
%components. Each graph is represented as an adjacency matrix.
for comp_i=1:length(comps)
    comp = comps{comp_i};
    nb_nodes = numel(comp);
    subgraph = zeros([nb_nodes, nb_nodes]);
    for i=1:numel(comp)
        nodeid = comp{i};
        
        
    end
end
end



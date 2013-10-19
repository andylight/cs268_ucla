function test_connected_components()
%TEST_CONNECTED_COMPONENTS Summary of this function goes here
%   Detailed explanation goes here
% G0 := Three connected components
G0 = [[1 1 1 0 0 0];
      [1 1 1 0 0 0];
      [1 1 1 0 0 0];
      [0 0 0 1 1 0];
      [0 0 0 1 1 0];
      [0 0 0 0 0 1]];

comps = connected_components(G0);
nb_comps = numel(comps);
disp(sprintf('nb_comps is: %d (should be: %d)', nb_comps, 3));

end


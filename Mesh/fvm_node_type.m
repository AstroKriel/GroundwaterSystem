classdef fvm_node_type < uint32
    % Author: Mitchell Johnson
    enumeration
        ALLUVIUM  (1),
        CONFINING (2),
        VOLCANICS (3),
        COAL      (4)
        % NEUTRAL (5),  % DO NOT NEED THESE
        % BOUNDARY (6), % DO NOT NEED THESE
    end
end
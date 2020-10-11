classdef fvm_node
    % Author: Mitchell Johnson
    properties
        node_type
    end
    
    methods
        function obj = fvm_node(node_type)
            obj.node_type = node_type;
        end
    end
end
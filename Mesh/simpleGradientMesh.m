function obj = simpleGradientMesh(mesh)
    % Author: Mitchell Johnson
    % Produces a simple non-uniform mesh to fit the unique material boundaries
    % and other constraints (eg, the river location) of our system.

    obj = struct;
    obj.base_width = mesh.c - 1;
    obj.base_height = mesh.r - 1;
    obj.physical_width = 500;
    obj.physical_height = 100;
    obj.x_delta = mesh.c_div; % 1
    obj.z_delta = mesh.r_div; % 1

    % Points of interest - hardcoded atm, will refactor into OO if time allows
    obj.x_weights = [0 5 30 35 50];
    obj.z_weights = [0 1 8 10 16 20]; % 16 is the river

    % Node width
    obj.n_width = (obj.base_width * obj.x_delta) + 1;
    obj.n_height = (obj.base_height * obj.z_delta) + 1;
    obj.n_count = obj.n_width * obj.n_height;

    obj.psi = ones(obj.n_width, obj.n_height);
    [x, z] = get_uniform_step_sizes(obj);
    [X, Z] = meshgrid(0:x:obj.base_width, 0:z:obj.base_height);
    X = X'; Z = Z';

    x_gradient = produce_gradient(obj.n_width, obj.x_weights);
    y_gradient = produce_gradient(obj.n_height, obj.z_weights);

    for z = 1:1:(obj.base_height * obj.z_delta) + 1
        X(:,z) = x_gradient;
    end
    for x = 1:1:(obj.base_width * obj.x_delta) + 1
        Z(x,:) = y_gradient;
    end

    % Save results
    obj.psi_x_loc = X;
    obj.psi_z_loc = Z;

    % Mesh
    obj.x_node_locations = obj.psi_x_loc(:,1)*(obj.physical_width / obj.psi_x_loc(end,1));
    obj.z_node_locations = flip(obj.psi_z_loc(1,:)*(obj.physical_height / obj.psi_z_loc(1,end)))';

    obj.deltax_e = zeros(obj.n_width, 1);
    for x = 1:1:obj.n_width-1
        obj.deltax_e(x) = (obj.x_node_locations(x+1) - obj.x_node_locations(x)) / 2;
    end
    obj.deltax_w = zeros(obj.n_width, 1);
    for x = 2:1:obj.n_width
        obj.deltax_w(x) = (obj.x_node_locations(x) - obj.x_node_locations(x-1)) / 2;
    end
    obj.deltaz_s = zeros(obj.n_height, 1);
    for x = 1:1:obj.n_height-1
        obj.deltaz_s(x) = (obj.z_node_locations(x) - obj.z_node_locations(x+1)) / 2;
    end
    obj.deltaz_n = zeros(obj.n_height, 1);
    for x = 2:1:obj.n_height
        obj.deltaz_n(x) = (obj.z_node_locations(x-1) - obj.z_node_locations(x)) / 2;
    end
    obj.Delta_x = zeros(obj.n_width, 1);
    for x = 1:1:obj.n_width
        obj.Delta_x(x) = obj.deltax_e(x) + obj.deltax_w(x);
    end
    obj.Delta_z = zeros(obj.n_height, 1);
    for x = 1:1:obj.n_height
        obj.Delta_z(x) = obj.deltaz_n(x) + obj.deltaz_s(x);
    end

    % Assign node types (hardcoded boundaries atm, can refactor to OO if time
    % allows)
    obj.m_split = zeros(obj.n_width * obj.n_height, 4);
    % obj.node_types = ones(obj.base_width * obj.x_delta, obj.base_height * obj.z_delta);
    for x = 1:1:obj.n_width
        for z = obj.n_height:-1:1
            ordering_no = ((x - 1) * obj.n_height) + (obj.n_height - z + 1);
            x_loc = obj.x_node_locations(x);
            z_loc = obj.z_node_locations(obj.n_height - z + 1);
            
            % 1 2
            % 3 4
            
            % Main Chunks
            if z_loc < 50
                obj.m_split(ordering_no, 1) = fvm_node_type.COAL;
                obj.m_split(ordering_no, 2) = fvm_node_type.COAL;
                obj.m_split(ordering_no, 3) = fvm_node_type.COAL;
                obj.m_split(ordering_no, 4) = fvm_node_type.COAL;
            end
            if x_loc < 300 && z_loc > 40
                obj.m_split(ordering_no, 1) = fvm_node_type.ALLUVIUM;
                obj.m_split(ordering_no, 2) = fvm_node_type.ALLUVIUM;
                obj.m_split(ordering_no, 3) = fvm_node_type.ALLUVIUM;
                obj.m_split(ordering_no, 4) = fvm_node_type.ALLUVIUM;
            end
            if x_loc > 300 && z_loc > 50
                obj.m_split(ordering_no, 1) = fvm_node_type.VOLCANICS;
                obj.m_split(ordering_no, 2) = fvm_node_type.VOLCANICS;
                obj.m_split(ordering_no, 3) = fvm_node_type.VOLCANICS;
                obj.m_split(ordering_no, 4) = fvm_node_type.VOLCANICS;
            end
            if x_loc > 50 && x_loc < 350 && z_loc > 40 && z_loc < 50
                obj.m_split(ordering_no, 1) = fvm_node_type.CONFINING;
                obj.m_split(ordering_no, 2) = fvm_node_type.CONFINING;
                obj.m_split(ordering_no, 3) = fvm_node_type.CONFINING;
                obj.m_split(ordering_no, 4) = fvm_node_type.CONFINING;
            end
            
            % Boundaries
            if x_loc > 0 && x_loc < 50 && z_loc == 40
                obj.m_split(ordering_no, 1) = fvm_node_type.ALLUVIUM;
                obj.m_split(ordering_no, 2) = fvm_node_type.ALLUVIUM;
                obj.m_split(ordering_no, 3) = fvm_node_type.COAL;
                obj.m_split(ordering_no, 4) = fvm_node_type.COAL;
            end
            if x_loc > 50 && x_loc < 350 && z_loc == 40
                obj.m_split(ordering_no, 1) = fvm_node_type.CONFINING;
                obj.m_split(ordering_no, 2) = fvm_node_type.CONFINING;
                obj.m_split(ordering_no, 3) = fvm_node_type.COAL;
                obj.m_split(ordering_no, 4) = fvm_node_type.COAL;
            end
            if x_loc == 350 && z_loc < 50 && z_loc > 40
                obj.m_split(ordering_no, 1) = fvm_node_type.CONFINING;
                obj.m_split(ordering_no, 2) = fvm_node_type.COAL;
                obj.m_split(ordering_no, 3) = fvm_node_type.CONFINING;
                obj.m_split(ordering_no, 4) = fvm_node_type.COAL;
            end
            if x_loc > 300 && x_loc < 350 && z_loc == 50
                obj.m_split(ordering_no, 1) = fvm_node_type.VOLCANICS;
                obj.m_split(ordering_no, 2) = fvm_node_type.VOLCANICS;
                obj.m_split(ordering_no, 3) = fvm_node_type.CONFINING;
                obj.m_split(ordering_no, 4) = fvm_node_type.CONFINING;
            end
            if x_loc > 50 && x_loc < 300 && z_loc == 50
                obj.m_split(ordering_no, 1) = fvm_node_type.ALLUVIUM;
                obj.m_split(ordering_no, 2) = fvm_node_type.ALLUVIUM;
                obj.m_split(ordering_no, 3) = fvm_node_type.CONFINING;
                obj.m_split(ordering_no, 4) = fvm_node_type.CONFINING;
            end
            if x_loc == 50 && z_loc < 50 && z_loc > 40
                obj.m_split(ordering_no, 1) = fvm_node_type.ALLUVIUM;
                obj.m_split(ordering_no, 2) = fvm_node_type.CONFINING;
                obj.m_split(ordering_no, 3) = fvm_node_type.ALLUVIUM;
                obj.m_split(ordering_no, 4) = fvm_node_type.CONFINING;
            end
            if x_loc > 350 && z_loc == 50
                obj.m_split(ordering_no, 1) = fvm_node_type.VOLCANICS;
                obj.m_split(ordering_no, 2) = fvm_node_type.VOLCANICS;
                obj.m_split(ordering_no, 3) = fvm_node_type.COAL;
                obj.m_split(ordering_no, 4) = fvm_node_type.COAL;
            end
            if x_loc == 300 && z_loc > 50
                obj.m_split(ordering_no, 1) = fvm_node_type.ALLUVIUM;
                obj.m_split(ordering_no, 2) = fvm_node_type.VOLCANICS;
                obj.m_split(ordering_no, 3) = fvm_node_type.ALLUVIUM;
                obj.m_split(ordering_no, 4) = fvm_node_type.VOLCANICS;
            end
            
            % Boundary Corners
            if x_loc == 0 && z_loc == 40
                obj.m_split(ordering_no, 1) = fvm_node_type.ALLUVIUM;
                obj.m_split(ordering_no, 2) = fvm_node_type.ALLUVIUM;
                obj.m_split(ordering_no, 3) = fvm_node_type.COAL;
                obj.m_split(ordering_no, 4) = fvm_node_type.COAL;
            end
            if x_loc == 50 && z_loc == 40
                obj.m_split(ordering_no, 1) = fvm_node_type.ALLUVIUM;
                obj.m_split(ordering_no, 2) = fvm_node_type.CONFINING;
                obj.m_split(ordering_no, 3) = fvm_node_type.COAL;
                obj.m_split(ordering_no, 4) = fvm_node_type.COAL;
            end
            if x_loc == 50 && z_loc == 50
                obj.m_split(ordering_no, 1) = fvm_node_type.ALLUVIUM;
                obj.m_split(ordering_no, 2) = fvm_node_type.ALLUVIUM;
                obj.m_split(ordering_no, 3) = fvm_node_type.ALLUVIUM;
                obj.m_split(ordering_no, 4) = fvm_node_type.CONFINING;
            end
            if x_loc == 300 && z_loc == 50
                obj.m_split(ordering_no, 1) = fvm_node_type.ALLUVIUM;
                obj.m_split(ordering_no, 2) = fvm_node_type.VOLCANICS;
                obj.m_split(ordering_no, 3) = fvm_node_type.CONFINING;
                obj.m_split(ordering_no, 4) = fvm_node_type.CONFINING;
            end
            if x_loc == 350 && z_loc == 50
                obj.m_split(ordering_no, 1) = fvm_node_type.VOLCANICS;
                obj.m_split(ordering_no, 2) = fvm_node_type.VOLCANICS;
                obj.m_split(ordering_no, 3) = fvm_node_type.CONFINING;
                obj.m_split(ordering_no, 4) = fvm_node_type.COAL;
            end
            if x_loc == 350 && z_loc == 40
                obj.m_split(ordering_no, 1) = fvm_node_type.CONFINING;
                obj.m_split(ordering_no, 2) = fvm_node_type.COAL;
                obj.m_split(ordering_no, 3) = fvm_node_type.COAL;
                obj.m_split(ordering_no, 4) = fvm_node_type.COAL;
            end
            if x_loc == 500 && z_loc == 50
                obj.m_split(ordering_no, 1) = fvm_node_type.VOLCANICS;
                obj.m_split(ordering_no, 2) = fvm_node_type.VOLCANICS;
                obj.m_split(ordering_no, 3) = fvm_node_type.COAL;
                obj.m_split(ordering_no, 4) = fvm_node_type.COAL;
            end
            if x_loc == 300 && z_loc == 100
                obj.m_split(ordering_no, 1) = fvm_node_type.ALLUVIUM;
                obj.m_split(ordering_no, 2) = fvm_node_type.VOLCANICS;
                obj.m_split(ordering_no, 3) = fvm_node_type.ALLUVIUM;
                obj.m_split(ordering_no, 4) = fvm_node_type.VOLCANICS;
            end
            
            % Domain Boundaries
            if x_loc == 0
                obj.m_split(ordering_no, 1) = NaN;
                obj.m_split(ordering_no, 3) = NaN;
            end
            if z_loc == 0
                obj.m_split(ordering_no, 3) = NaN;
                obj.m_split(ordering_no, 4) = NaN;
            end
            if x_loc == obj.physical_width
                obj.m_split(ordering_no, 2) = NaN;
                obj.m_split(ordering_no, 4) = NaN;
            end
            if z_loc == obj.physical_height
                obj.m_split(ordering_no, 1) = NaN;
                obj.m_split(ordering_no, 2) = NaN;
            end
        end
    end

end

function [x, y] = get_uniform_step_sizes(obj)
    % TODO: simplify
    x = (obj.base_width/((obj.base_width * obj.x_delta)));
    y = (obj.base_height/((obj.base_height * obj.z_delta)));
end

function gradient_x = produce_gradient(count, sinks)
    count = count - 1;
    [~, sink_size] = size(sinks);
    placeable = count;
    % Get aim-to-fit-between nodes
    % If we have 5 nodes, we have 4 gaps, for example.
    number_want_to_fit_in_gaps = floor(placeable / (sink_size - 1));
    number_want_to_fit_in_gaps_leftovers = mod(placeable, (sink_size - 1));
    % List available to place by gap
    available_list = ones(1, sink_size - 1) * number_want_to_fit_in_gaps;
    % Fill in the biggest gaps with leftovers
    [~, I] = maxk(diff(sinks), number_want_to_fit_in_gaps_leftovers);
    for j = 1:size(I, 2)
        available_list(I(j)) = available_list(I(j)) + 1;
    end
    % Produce initial map
    gradient_x = [];
    for j = 1:size(sinks,2)-1
        gradient_x = [gradient_x produce_quadratic(sinks(j), sinks(j+1), available_list(j))];
    end
    gradient_x = [gradient_x sinks(end)];
end

function qx_loc = produce_quadratic(start, finish, count)
    quadratic = @(x) -4*x.^2 + 4*x;
    x_loc = (-(finish-start)/2):((finish-start)/count):((finish-((finish-start)/count))-start-((finish-start)/2));
    if (size(x_loc, 2) ~= count)
        x_loc = [x_loc ((finish-((finish-start)/count))-start-((finish-start)/2))];
    end
    % Display size for debugging
    % size(x_loc, 2)
    x_norm = quadratic(0:1/count:1-(1/count));
    x_loc = x_loc + (0.25 .* x_loc .* x_norm);
    qx_loc = x_loc + (finish-start)/2 + start;
end
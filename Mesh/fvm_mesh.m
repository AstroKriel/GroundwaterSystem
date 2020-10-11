%% MESH class

% Author: Mitchell Johnson
% Purpose: Provide an object-oriented class to handle the mesh
% discretization within the context of our groundwater modelling problem.
% This was originally designed to handle gradient mesh as well, however
% the non-uniform functionality was offloaded to simpleGradientMesh to
% make it easier to understand.

% Features:
% - Generate a uniform mesh with provided x, z, relative step distance
% - Generate a uniform mesh with provided x, z, and 0-100% randomness
% - Graph itself (in whatever node points are provided)
% - Node points need to have x, z, type of material, and other constants

%% Class Implementation
classdef fvm_mesh
    % Author: Mitchell Johnson
    
    % ===============
    % Constants
    % ===============
    % These variables are any constants used in our discretization. 
    properties (Constant)
        DISPLAY_LOGGING = 0 % 0 = no, 1 = some, 2 = all
        NODE_GRID_COLOR = [0 0 0] + 0.05*18
    end
    
    % ===============
    % Private Properties
    % ===============
    % These are state variable that should not be accessible outside of the
    % fvm_mesh class
    properties (SetAccess = private)
        % Our current position in time
        time_step = 0
        
        % Physical Size
        width
        height
        
        % Fractional step-size of mesh in x-direction
        x_delta = 1
        % Fractional step-size of mesh in z-direction
        z_delta = 1
        
        % Mesh type
        mesh_type
        
        % COMPLETED STATE
        IS_BAKED = 0
        
        % The mesh (x, z) locations
        psi_x_loc
        psi_z_loc
                
        % The material types around each node
        node_types = [];
        
        % Weights for computing the mesh non-uniformity
        x_weights = [];
        y_weights = [];
    end
    
    % ===============
    % Variable Properties
    % ===============
    % These variables contain the core state of our discretization. They
    % are publically modifiable and accessible
    properties
        % Our maximum time (beyond which time_step cannot proceed)
        maximum_time = 50* 365 * 24 * 60 * 60
        % The default time (in seconds) to jump forward for every step.
        time_step_delta = 60
        % The percentage of randomness in the x-direction
        x_randomness = 0
        % The percentage of randomness in the z-direction
        z_randomness = 0
        
        % OPTIONS
        SHOW_NODES = 1
        SHOW_NODE_GRID = 0
        SHOW_VOLUME_MESH = 1
        
        % psi values (the actual discretization)
        psi
        
        % deltax_e is a vector that is Nx1 that gives the distance from the
        % current node P to the east control volume face.
        deltax_e
        % deltax_w is a vector that is Nx1 that gives the distance from the
        % current node P to the west control volume face.
        deltax_w
        % deltaz_n is a vector that is Nx1 that gives the distance from the
        % current node P to the north control volume face.
        deltaz_n
        % deltaz_s is a vector that is Nx1 that gives the distance from the
        % current node P to the south control volume face.
        deltaz_s
        
        % Delta_x is the distance between the nodes. If we want the
        % distance between i and i+1, take Delta_x(i). If instead we want
        % the distance between i-1 and i take Delta_x(i-1). 
        Delta_x
        
        % Delta_z is the distance between the nodes. If we want the
        % distance between i and the i+1, take Delta_z(i). If instead we
        % want the distance between i-1 and i take Delta_z(i-1).
        Delta_z
        
        % x_node_locations shows the x value of the column of nodes. Is an
        % r x 1 vector.
        x_node_locations
        % z_node_locations shows the z value of the row of nodes. Is a c x
        % 1 vector. 
        z_node_locations
        
        physical_width
        physical_height
    end
    
    % ===============
    % Dependent Properties
    % ===============
    % These variables depend at runtime upon the values of the other
    % state variables. We define them here so that they are re-computed at
    % runtime whenever they are accessed.
    properties (Dependent)
        % The number of nodes along the x-direction of the mesh
        n_width
        % The number of nodes along the z-direction of the mesh
        n_height
        % The total number of nodes
        n_count
        
        % material types around each node
        m_split
    end
    
    % ===============
    % Methods
    % ===============
    % These functions allow us to modify the state of our discretization
    % within the boundaries of defined behaviour, as well as run stateless
    % behaviours such as displaying plots and outputting logging to the
    % Command Window.
    methods
        % Constructor
        function obj = fvm_mesh(mesh_type, width, height, x_delta, z_delta, p_width, p_height)
            if nargin == 3
                obj.mesh_type = mesh_type;
                obj.width = width;
                obj.height = height;
                obj.x_delta = 1;
                obj.z_delta = 1;
                obj.node_types = ones(width * x_delta, height * z_delta);
                obj.physical_width = width; 
                obj.physical_height = height;
            elseif nargin == 7
                obj.mesh_type = mesh_type;
                obj.width = width;
                obj.height = height;
                obj.x_delta = x_delta;
                obj.z_delta = z_delta;
                obj.node_types = ones(width * x_delta, height * z_delta);
                obj.physical_width = p_width;
                obj.physical_height = p_height;
            else
                error('Invalid # of arguments. fvm_mesh must take 3 or 7');
            end
        end
        
        % Get Methods
        % (these compute the values of the dependent properties at runtime)
        function val = get.n_width(obj)
            if obj.IS_BAKED
                [val, ~] = size(obj.psi);
            else
                error('Must finalize mesh before requesting width');
            end
        end
        function val = get.n_height(obj)
            if obj.IS_BAKED
                [~, val] = size(obj.psi);
            else
                error('Must finalize mesh before requesting height');
            end
        end
        function val = get.n_count(obj)
           val = obj.n_width * obj.n_height; 
        end
        
        function val = get.x_node_locations(obj)
            val = obj.psi_x_loc(:,1)*(obj.physical_width / obj.psi_x_loc(end,1));
        end
        function val = get.z_node_locations(obj)
            val = flip(obj.psi_z_loc(1,:)*(obj.physical_height / obj.psi_z_loc(1,end)))';
        end
        function val = get.deltax_e(obj)
            val = zeros(obj.n_width, 1);
            for x = 1:1:obj.n_width-1
                val(x) = (obj.x_node_locations(x+1) - obj.x_node_locations(x)) / 2;
            end
        end
        function val = get.deltax_w(obj)
            val = zeros(obj.n_width, 1);
            for x = 2:1:obj.n_width
                val(x) = (obj.x_node_locations(x) - obj.x_node_locations(x-1)) / 2;
            end
        end
        function val = get.deltaz_s(obj)
            val = zeros(obj.n_height, 1);
            for x = 1:1:obj.n_height-1
                val(x) = (obj.z_node_locations(x) - obj.z_node_locations(x+1)) / 2;
            end
        end
        function val = get.deltaz_n(obj)
            val = zeros(obj.n_height, 1);
            for x = 2:1:obj.n_height
                val(x) = (obj.z_node_locations(x-1) - obj.z_node_locations(x)) / 2;
            end
        end
        function val = get.Delta_x(obj)
            val = zeros(obj.n_width, 1);
            for x = 1:1:obj.n_width
                val(x) = obj.deltax_e(x) + obj.deltax_w(x);
            end
        end
        function val = get.Delta_z(obj)
            val = zeros(obj.n_height, 1);
            for x = 1:1:obj.n_height
                val(x) = obj.deltaz_n(x) + obj.deltaz_s(x);
            end
        end
        
        function val = get.m_split(obj)
            % Define an empty array with the correct shape
            val = ones(obj.n_width * obj.n_height, 4);
            
            % Iterate through and fill with data appropriately
            % (using a hard-coded node-ordering atm)
            for x = 1:1:obj.n_width
                for y = obj.n_height:-1:1
                    ordering_no = ((x - 1) * obj.n_height) + (obj.n_height - y + 1);
                    
                    if x > 1 && y < obj.n_height
                        val(ordering_no, 1) = obj.node_types(x - 1, y);
                    else
                        val(ordering_no, 1) = NaN;
                    end
                    
                    if x < obj.n_width && y < obj.n_height
                        val(ordering_no, 2) = obj.node_types(x, y);
                    else
                        val(ordering_no, 2) = NaN;
                    end
                    
                    if x > 1 && y > 1
                        val(ordering_no, 3) = obj.node_types(x - 1, y - 1);
                    else
                        val(ordering_no, 3) = NaN;
                    end
                    
                    if x < obj.n_width && y > 1
                        val(ordering_no, 4) = obj.node_types(x, y - 1);
                    else
                        val(ordering_no, 4) = NaN;
                    end
                    
                    % Uncomment to add debug columns
                    % val(ordering_no, 5) = x;
                    % val(ordering_no, 6) = y;
                end
            end
            
            % val = ones(2, 2);
        end
        
        % Set Methods
        % (We use these to ensure type consistency)
        
        
        % Computation Functions
        function obj = step(obj, amount_of_time)
            % TODO: implement this with PDE and temporal discretization
        end
        
        % Utility Functions
        function obj = define_region(obj, x_start, z_start, x_end, z_end, region_type)
            % Scale to fineness
            x_start = (x_start * obj.x_delta);
            z_start = (z_start * obj.z_delta);
            x_end = (x_end * obj.x_delta);
            z_end = (z_end * obj.z_delta);
            % Modify node types matrix
            obj.node_types([x_start+1:x_end],[z_start+1:z_end]) = region_type;
            % Update boundary vectors
            if ~ismember(x_start, obj.x_weights)
                obj.x_weights = sort([obj.x_weights, x_start]);
            end
            if ~ismember(x_end, obj.x_weights)
                obj.x_weights = sort([obj.x_weights, x_end]);
            end
            if ~ismember(z_start, obj.y_weights)
                obj.y_weights = sort([obj.y_weights, z_start]);
            end
            if ~ismember(z_end, obj.y_weights)
                obj.y_weights = sort([obj.y_weights, z_end]);
            end
        end
        
        function obj = bake(obj)
            % Set the mesh to a matrix of ones, with correct step size.
            obj.psi = ones((obj.width * obj.x_delta) + 1, (obj.height * obj.z_delta) + 1);

            % Produce a uniform meshgrid for the actual locations
            [x, z] = obj.get_uniform_step_sizes(obj);
            [X, Z] = meshgrid(0:x:obj.width, 0:z:obj.height);
            X = X'; Z = Z';
            
            % Gradient? If so, override the uniform.
            % TODO: repair for new fineness
            if obj.mesh_type == fvm_mesh_type.GRADIENT
                x_gradient = obj.produce_gradient(obj, (obj.width * obj.x_delta)+1, obj.x_weights);
                y_gradient = obj.produce_gradient(obj, (obj.height * obj.z_delta)+1, obj.y_weights);
                
                for z = 1:1:(obj.height * obj.z_delta) + 1
                    X(:,z) = x_gradient;
                end
                for x = 1:1:(obj.width * obj.x_delta) + 1
                    Z(x,:) = y_gradient;
                end
            end
            
            % Add noise?
%             if obj.mesh_type == fvm_mesh_type.UNIFORM_with_granular_noise
%                 x_noise = 0.15*x*randn(size(X));
%                 z_noise = 0.15*z*randn(size(Z));
%                 for j = 1:1:size(obj.x_weights,2)
%                     x_noise(obj.x_weights(j)+1,:) = 0;
%                     z_noise(obj.x_weights(j)+1,:) = 0;
%                 end
%                 for j = 1:1:size(obj.y_weights,2)
%                     z_noise(:,obj.y_weights(j)+1) = 0;
%                     x_noise(:,obj.y_weights(j)+1) = 0;
%                 end
%                 X = X + x_noise;
%                 Z = Z + z_noise;
%             end

            obj.psi_x_loc = X;
            obj.psi_z_loc = Z;
            
            % Bake the mesh
            obj.IS_BAKED = 1;
        end
        
        function plot(obj)
            % Setup figure and force pretty axis size
            figure('Renderer', 'painters', 'Position', [10 40 1700 1700*(20/50)])
            axis([0 obj.width 0 obj.height])
            hold on;
            
            im = imagesc([0.5 / obj.x_delta obj.width - (0.5 / obj.x_delta)], [0.5 / obj.z_delta obj.height - (0.5 / obj.z_delta)], obj.node_types');
            im.AlphaData = .2;
            
            % Loop through our points
            for x = 1:1:obj.n_width
                for y = 1:1:obj.n_height
                    xn = obj.psi_x_loc(x, y);
                    yn = obj.psi_z_loc(x, y);
                    if obj.SHOW_NODE_GRID
                        if x ~= obj.n_width
                            plot([xn;obj.psi_x_loc(x+1, y)],[yn;obj.psi_z_loc(x+1, y)],'linewidth',1,'Color',obj.NODE_GRID_COLOR);
                        end
                        if y ~= obj.n_height
                            plot([xn;obj.psi_x_loc(x, y+1)],[yn;obj.psi_z_loc(x, y+1)],'linewidth',1,'Color',obj.NODE_GRID_COLOR);
                        end
                    end
                    if obj.SHOW_NODES
                        plot(xn, yn, 'b.');
                    end
                end
            end
            
            % Volume mesh logic
            volume_mesh_x = ones(obj.n_width - 1, obj.n_height - 1);
            volume_mesh_y = ones(obj.n_width - 1, obj.n_height - 1);
            [ux, uy] = obj.get_uniform_step_sizes(obj);
            if obj.SHOW_VOLUME_MESH
                % Pre-compute mesh
                for x = 1:1:obj.n_width-1
                    for y = 1:1:obj.n_height-1
                        %    t--s
                        %    |  |
                        % l--n--r
                        %    |
                        %    b
                        xn = obj.psi_x_loc(x, y);
                        yn = obj.psi_z_loc(x, y);
                        % Get Important Points Around
                        xr = obj.psi_x_loc(x+1, y);
                        yr = obj.psi_z_loc(x+1, y);
                        xt = obj.psi_x_loc(x, y+1);
                        yt = obj.psi_z_loc(x, y+1);
                        xs = obj.psi_x_loc(x+1, y+1);
                        ys = obj.psi_z_loc(x+1, y+1);
                        % Compute volume mesh point
                        [vm, vy] = obj.intersect_four(xn, yn, xs, ys, xt, yt, xr, yr);
                        volume_mesh_x(x, y) = vm;
                        volume_mesh_y(x, y) = vy;
                    end
                end
                % Display mesh
                for x = 1:1:obj.n_width-1
                    for y = 1:1:obj.n_height-1
                        xn = volume_mesh_x(x, y);
                        yn = volume_mesh_y(x, y);
                        if x ~= obj.n_width-1
                            plot([xn;volume_mesh_x(x+1, y)],[yn;volume_mesh_y(x+1, y)],'linewidth',1,'Color',[0.8 0 0.9]);
                            if x == 1
                                y_b = (obj.psi_z_loc(x,y) + obj.psi_z_loc(x,y+1)) / 2;
                                plot([0;xn],[y_b;yn],'linewidth',1,'Color',[0.8 0 0.9]);
                            end
                        else
                            if x == obj.n_width-1
                                y_b = (obj.psi_z_loc(x,y) + obj.psi_z_loc(x,y+1)) / 2;
                                plot([xn;x*ux],[yn;y_b],'linewidth',1,'Color',[0.8 0 0.9]);
                            end
                        end
                        if y ~= obj.n_height-1
                            plot([xn;volume_mesh_x(x, y+1)],[yn;volume_mesh_y(x, y+1)],'linewidth',1,'Color',[0.8 0 0.9]);
                            if y == 1
                                x_b = (obj.psi_x_loc(x,y) + obj.psi_x_loc(x+1,y)) / 2;
                                plot([x_b;xn],[0;yn],'linewidth',1,'Color',[0.8 0 0.9]);
                            end
                        else
                            if y == obj.n_height-1
                                x_b = (obj.psi_x_loc(x,y) + obj.psi_x_loc(x+1,y)) / 2;
                                plot([xn;x_b],[yn;y*uy],'linewidth',1,'Color',[0.8 0 0.9]);
                            end
                        end
                    end
                end
            end
            
            hold off;
        end
    end
    
    % ===============
    % Static Methods
    % ===============
    % These methods exist outside the context of the class and allow for
    % stateless objects to be constructed, etc. They do not require state,
    % and only exist to act as utility functions.
    methods (Static)
        
        function [x, y] = get_uniform_step_sizes(obj)
            % TODO: can this be simple division?
            x = (obj.width/((obj.width * obj.x_delta)));
            y = (obj.height/((obj.height * obj.z_delta)));
        end
        
        function [x, y] = intersect_four(x1, y1, x2, y2, x3, y3, x4, y4)
            % Computes the midpoint intersection of all points
            a1 = y2 - y1;
            b1 = x1 - x2;
            c1 = a1 * x1 + b1 * y1;
            a2 = y4 - y3;
            b2 = x3 - x4;
            c2 = a2 * x3 + b2 * y3;
            det = a1 * b2 - a2 * b1;
            x = (b2 * c1 - b1 * c2) / det; 
            y = (a1 * c2 - a2 * c1) / det; 
        end
        
        function gradient_x = produce_gradient(obj, count, sinks)
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
                gradient_x = [gradient_x obj.produce_quadratic(sinks(j), sinks(j+1), available_list(j))];
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
        
    end
    
end
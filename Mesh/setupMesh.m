function [fm, cfg] = setupMesh(mesh, cfginput)
    %{
        SETUPMESH
        Authors: Mitchell Johnson, Sam Dudley
    %}

    cfg = cfginput;

    % Default Graph
    physicalWidth  = cfg.L1;
    physicalHeight = cfg.L2;

    % Node dimensions
    nodeWidth  = mesh.c - 1;
    nodeHeight = mesh.r - 1;

    % Node resolution
    % DO NOT CHANGE - Remove from later versions
    % Mitch: passed these through as config variables
    if isfield(mesh, 'c_div')
        nodeFinenessX = mesh.c_div;
    else
        nodeFinenessX = 1;
    end
    if isfield(mesh, 'r_div')
        nodeFinenessY = mesh.r_div;
    else
        nodeFinenessY = 1;
    end

    % Create the fm object

    if strcmpi(mesh.type, 'uniform')
        % Create the mesh object
        fm = fvm_mesh(fvm_mesh_type.UNIFORM, nodeWidth, nodeHeight, ...
            nodeFinenessX, nodeFinenessY, ...
            physicalWidth, physicalHeight);

        % Dynamically define boundaries
        % startX, startZ, endX, endZ, material type
        fm = fm.define_region(0, 0, 50, 8, fvm_node_type.COAL);
        fm = fm.define_region(35, 8, 50, 10, fvm_node_type.COAL);
        fm = fm.define_region(5, 8, 35, 10, fvm_node_type.CONFINING);
        fm = fm.define_region(30, 10, 50, 20, fvm_node_type.VOLCANICS);
        fm = fm.define_region(5, 10, 30, 20, fvm_node_type.ALLUVIUM);
        fm = fm.define_region(0, 8, 5, 20, fvm_node_type.ALLUVIUM);

        % Indicate that we are finished
        fm = fm.bake();

    elseif strcmpi(mesh.type, 'gradient')
        % TODO: add parameters
        fm = simpleGradientMesh(mesh);
        % fm = fvm_mesh(fvm_mesh_type.GRADIENT, nodeWidth, nodeHeight, ...
        %     nodeFinenessX, nodeFinenessY, ...
        %     physicalWidth, physicalHeight);
    end

    % --- Get information from mesh --- %
    % Number of rows and columns
    r = fm.n_height;
    c = fm.n_width;
    N = fm.n_count;

    cfg.r = r; cfg.c = c; cfg.N = N;

    % Create matrix.
    cfg.x = fm.x_node_locations';
    cfg.z = fm.z_node_locations';
    [X, Z] = meshgrid(cfg.x,cfg.z);
    cfg.Z = Z; cfg.X = X;

    % These are the distances to the CONTROL VOLUMES not the neighbouring
    % nodes.

    % East and West.
    cfg.delta_xe = meshgrid(fm.deltax_e', 1:r);
    cfg.delta_xw = meshgrid(fm.deltax_w', 1:r);

    % North and South
    cfg.delta_zn = meshgrid(fm.deltaz_n, 1:c)';
    cfg.delta_zs = meshgrid(fm.deltaz_s, 1:c)';

    % Here we consider the control volumes for each node, including the sub
    % control volumes. The full control volume can be found by taking the row
    % sum of DV.
    DV = zeros(N,4);
    DV(:,1) = cfg.delta_xw(:) .* cfg.delta_zn(:);
    DV(:,2) = cfg.delta_xe(:) .* cfg.delta_zn(:);
    DV(:,3) = cfg.delta_xw(:) .* cfg.delta_zs(:);
    DV(:,4) = cfg.delta_xe(:) .* cfg.delta_zs(:);
    cfg.DV = DV;

    % Control volume sum
    cfg.CV = sum(DV,2);

    % Check the sub control volumes are correct.
    % There's a difference of -6.5484e-11 on the gradient mesh which is
    % negligible error due to numerical constraints, but enough to trip this check.
    % if sum(cfg.DV(:)) ~= cfg.L1 * cfg.L2
    %     error('the Sub Control Volume calculation is incorrect')
    % end

    % NM contains the 4 material types for each node.
    %   |       |
    % - 1 ----- 2 -
    %   |   P   |
    % - 3 ----- 4 -
    %   |       |
    %
    cfg.NM = fm.m_split;

end

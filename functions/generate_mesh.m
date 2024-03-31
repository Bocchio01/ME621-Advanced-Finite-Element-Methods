function [coordinates, connectivity] = generate_mesh(data, fem)

Lx_element = data.Lx / fem.Nx;
Ly_element = data.Ly / fem.Ny;

N_nodes = (fem.Nx + 1) * (fem.Ny + 1);
N_elements = fem.Nx * fem.Ny;

coordinates = zeros(N_nodes, fem.N_DOF_node);
connectivity = zeros(N_elements, fem.N_nodes_element);

el_idx = 0;
for i = 1:(fem.Nx + 1)
    for j = 1:(fem.Ny + 1)
        el_idx = el_idx + 1;
        coordinates(el_idx, 1) = Lx_element * (i - 1);
        coordinates(el_idx, 2) = Ly_element * (j - 1);
    end
end

el_idx = 0;
for i = 1:fem.Nx
    for j = 1:fem.Ny
        el_idx = el_idx + 1;
        connectivity(el_idx, 1) = el_idx + floor((el_idx-1)/fem.Ny);
        connectivity(el_idx, 2) = el_idx + floor((el_idx-1)/fem.Ny) + fem.Ny + 1;
        connectivity(el_idx, 3) = el_idx + floor((el_idx-1)/fem.Ny) + fem.Ny + 2;
        connectivity(el_idx, 4) = el_idx + 1 + floor((el_idx-1)/fem.Ny);
    end
end

end


function L = compute_gather_matrix(element_idx, N_nodes_element, N_nodes)

condition = @(i, j) ((N_nodes_element - 1) * (element_idx - 1) + i == j);

L = zeros(N_nodes_element, N_nodes);

for i = 1:N_nodes_element
    for j = 1:N_nodes
        L(i, j) = 1 * condition(i, j) + 0 * ~condition(i, j);
    end
end


end
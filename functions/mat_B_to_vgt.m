function vgt_B = mat_B_to_vgt(mat_B, N_nodes_element, N_DOF_node)

vgt_B = zeros(3, N_nodes_element * N_DOF_node);

for node_idx = 1:N_nodes_element
    vgt_B(:, (N_DOF_node*node_idx-1):(N_DOF_node*node_idx)) = [mat_B(1, node_idx), 0; 0, mat_B(2, node_idx); mat_B(2, node_idx), mat_B(1, node_idx)];
end

end
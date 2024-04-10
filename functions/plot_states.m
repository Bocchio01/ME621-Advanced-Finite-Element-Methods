function figure_results_analysis = plot_states(mesh, states, t)

figure_results_analysis = figure('Name', 'Results analysis', 'NumberTitle', 'off');
tiles = tiledlayout(1, 3);

if nargin == 3
    title(tiles, ['t = ' num2str(t) ' s']);
end

% Initial vs. Deformed configuration
nexttile([1, 2]);
hold on
grid on
axis equal;

for element_idx = 1:size(mesh.connectivity, 1)
    
    connections = mesh.connectivity(element_idx, [1:end, 1]);
    idx_nodes = reshape((2 * connections.' + [-1, 0])', 1, []);
    
    u_mesh = mesh.coordinates(connections, :);
    u_deformed = states.global.u(idx_nodes, end);
    
    plot(u_mesh(:, 1), u_mesh(:, 2), '-ob', 'LineWidth', 2)
    plot(u_deformed(1:2:end), u_deformed(2:2:end), '-or', 'LineWidth', 2)
    
end

title('Initial vs. Deformed configuration');
legend('Initial configuration', 'Deformed configuration', 'Location', 'Best')
xlabel('x [m]');
ylabel('y [m]');


% Shear stress vs. Shear strain
nexttile(3);
hold on
grid on

plot(states.targets.gamma, states.targets.sigma12 * 1e-6, '-', 'LineWidth', 2)

title('Shear stress vs. Shear strain');
xlabel('Shear strain [rad]');
ylabel('Shear stress [MPa]');

end

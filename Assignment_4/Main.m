% ME621 - Advanced Finite Element Methods
% Assignment 4
% Tommaso Bocchietti
% A.Y. 2023/24 - W24

% vgt_ prefix indicate that the variable is written in Voigt notation
% mat_ prefix indicate that the variable is written in Matrix notation

clc
clear variables
close all


%% Problem datas

data = struct( ...
    'Lx', 1, ...
    'Ly', 1, ...
    'rho', 2700, ...
    'E', 70e+9, ...
    'nu', 0.3, ...
    'vx_top', 1);

fem = struct( ...
    'N_nodes_element', 4, ...
    'N_DOF_node', 2, ...
    'N_DOF', NaN, ...
    'Nx', 1, ...
    'Ny', 1, ...
    'constitutive_model', "Trusdell");

output = struct('plots', true);

%% Mesh generation

mesh = struct();
[mesh.coordinates, mesh.connectivity] = generate_mesh(data, fem);
fem.N_DOF = numel(mesh.coordinates);


%% Boundary condition setup

idx_nodes_bottom = reshape((2 * find(mesh.coordinates(:, 2) == 0) + [-1, 0])', 1, []);
idx_nodes_top = reshape((2 * find(mesh.coordinates(:, 2) == data.Ly) + [-1, 0])', 1, []);

initial_coordinates= reshape(mesh.coordinates.', 1, [])';
boundary = struct( ...
    'displacement', struct( ...
    'index', [idx_nodes_bottom idx_nodes_top(2:2:end)], ...
    'value', initial_coordinates([idx_nodes_bottom idx_nodes_top(2:2:end)]) ...
    ), ...
    'velocities', struct( ...
    'index', [idx_nodes_bottom idx_nodes_top(2:2:end) idx_nodes_top(1:2:end)], ...
    'value', [zeros(1, length([idx_nodes_bottom idx_nodes_top(2:2:end)])) data.vx_top * ones(1, length(idx_nodes_top(1:2:end)))] ...
    ));

clear idx_nodes_bottom idx_nodes_top initial_coordinates


%% Solver structure definitions

solver = struct( ...
    'dt', zeros(1), ...
    'global', struct(), ...
    'integration', struct(), ...
    'element', struct(), ...
    'mat', struct(), ...
    'vgt', struct());

% solver.global
% Here we just initialize the vectors of position and velocities
% Boundary condition are not applied yet
solver.global.u(:, 1)      = reshape(mesh.coordinates.', 1, [])';
solver.global.v_half(:, 1) = zeros(fem.N_DOF, 1);

% solver.integration
solver.integration.weights = [5/9, 8/9, 5/9];
solver.integration.roots   = [-sqrt(3/5), 0, sqrt(3/5)];
% solver.integration.weights = [2];
% solver.integration.roots   = [0];

% solver.element
% Other fields of this struct will be populated directly inside the main loop cycle
solver.element.stress = zeros( ...
    size(mesh.connectivity, 1), ...
    length(solver.integration.roots), ...
    length(solver.integration.roots), ...
    3);

% solver.mat
solver.mat.N      = @(Xi, Eta)    1/4 * [(1-Xi)*(1-Eta), (1+Xi)*(1-Eta), (1+Xi)*(1+Eta), (1-Xi)*(1+Eta)];
solver.mat.dN     = @(Xi, Eta)    1/4 * [-1+Eta, 1-Eta, 1+Eta, -1-Eta; -1+Xi, -1-Xi, 1+Xi, 1-Xi];
solver.mat.F      = @(Xi, Eta, u) u * solver.mat.dN(Xi, Eta)';
solver.mat.J      = @(Xi, Eta, u) det(solver.mat.F(Xi, Eta, u));
solver.mat.B      = @(Xi, Eta, u) inv(solver.mat.F(Xi, Eta, u))' * solver.mat.dN(Xi, Eta);
solver.mat.stress = @(vgt_stress) vgt_stress_to_mat(vgt_stress);
solver.mat.strain = @(vgt_strain) vgt_strain_to_mat(vgt_strain);

% solver.vgt
solver.vgt.N      = @(N)          mat_N_to_vgt(N, fem.N_DOF_node);
solver.vgt.B      = @(B)          mat_B_to_vgt(B, fem.N_nodes_element, fem.N_DOF_node);
solver.vgt.stress = @(mat_stress) mat_stress_to_vgt(mat_stress);
solver.vgt.strain = @(mat_strain) mat_strain_to_vgt(mat_strain);
solver.vgt.C      = @(E, nu)      E / (1-nu^2) * [1, nu, 0; nu, 1, 0; 0, 0, (1-nu)/2]; % Plane stress model

solver.vgt.M = zeros(fem.N_DOF, fem.N_DOF);
M = zeros(fem.N_DOF_node * fem.N_nodes_element);

% Here element are the same in initial configuration and we can compute one single element mass matrix knowing
% that it will be the same for all the elements
for i = 1:length(solver.integration.roots)
    for j = 1:length(solver.integration.roots)
        
        % Variables for the loop
        idx_nodes = reshape((2 * mesh.connectivity(1, :).' + [-1, 0])', 1, []);
        u = reshape(solver.global.u(idx_nodes, 1), [fem.N_DOF_node, fem.N_nodes_element]);
        
        Xi = solver.integration.roots(i);
        Eta = solver.integration.roots(j);
        weight_Xi = solver.integration.weights(i);
        weight_Eta = solver.integration.weights(j);
        
        N = solver.mat.N(Xi, Eta);
        J = solver.mat.J(Xi, Eta, u);
        
        % Loop body
        M = M + weight_Xi * weight_Eta * data.rho * solver.vgt.N(N)' * solver.vgt.N(N) * J;
        
    end
end

for element_idx = 1:size(mesh.connectivity, 1)
    
    idx_nodes = reshape((2 * mesh.connectivity(element_idx, :).' + [-1, 0])', 1, []);
    solver.vgt.M(idx_nodes, idx_nodes) = solver.vgt.M(idx_nodes, idx_nodes) + M;
    
end

solver.vgt.M = diag(sum(solver.vgt.M, 2));

clear i j element_idx
clear idx_nodes u Xi Eta weight_Xi weight_Eta N J M


%% Plasticity model

plasticity = struct();

plasticity.Sy           = @(S)               sqrt(3/2 * sum(S .* S, 'all'));
% plasticity.hardening    = @(strain_p_bar)    200e6 + 325e6 * (strain_p_bar);
plasticity.hardening    = @(strain_p_bar)    200e6 + 325e6 * (strain_p_bar)^0.125;
plasticity.strain_p_bar = @(mat_strain_p)    sqrt(2/3 * sum(mat_strain_p .* mat_strain_p, 'all'));
plasticity.mat_strain_p = @(S, strain_p_bar) 3/2 * strain_p_bar / plasticity.Sy(S) * S;
plasticity.n_direction  = @(S)               3/2 * S / plasticity.Sy(S);

plasticity.phi          = @(S, strain_p_bar) plasticity.Sy(S) - plasticity.hardening(strain_p_bar);

%% Main loop

plotting_struct = struct( ...
    'sigma12', [0], ...
    'gamma', [0]);

n = 1;
t = 0;
C = solver.vgt.C(data.E, data.nu);

gamma(n) = [0];
solver.dt(n) = 5e-4;

tic
while (gamma(n) < 1)
    
    solver.global.f(:, n) = zeros(fem.N_DOF, 1);
    dt_critical = +inf;
    
    % 05) Stress update algorithm
    for element_idx = 1:size(mesh.connectivity, 1)
        
        idx_nodes = reshape((2 * mesh.connectivity(element_idx, :).' + [-1, 0])', 1, []);
        
        solver.element.vgt.f_int = zeros(fem.N_DOF_node * fem.N_nodes_element, 1);
        solver.element.vgt.f_ext = zeros(fem.N_DOF_node * fem.N_nodes_element, 1);
        
        solver.element.vgt.u = solver.global.u(idx_nodes, max(n-1, 1));
        solver.element.vgt.v_half = solver.global.v_half(idx_nodes, max(n-1, 1));
        
        solver.element.mat.u = reshape(solver.element.vgt.u, [fem.N_DOF_node, fem.N_nodes_element]);
        solver.element.mat.v_half = reshape(solver.element.vgt.v_half, [fem.N_DOF_node, fem.N_nodes_element]);
        
        if (n > 1)
            for i = 1:length(solver.integration.roots)
                for j = 1:length(solver.integration.roots)
                    
                    % Variables for the loop
                    Xi = solver.integration.roots(i);
                    Eta = solver.integration.roots(j);
                    weight_Xi = solver.integration.weights(i);
                    weight_Eta = solver.integration.weights(j);
                    
                    N = solver.mat.N(Xi, Eta);
                    B = solver.mat.B(Xi, Eta, solver.element.mat.u);
                    J = solver.mat.J(Xi, Eta, solver.element.mat.u);
                    
                    vgt_stress = squeeze(solver.element.stress(element_idx, i, j, :));
                    mat_stress = solver.mat.stress(vgt_stress);
                    
                    % Elastic predictor
                    switch fem.constitutive_model
                        case "Cauchy"
                            gamma_dot = (plotting_struct.gamma(n-1) - plotting_struct.gamma(max(n-2, 1))) / solver.dt(n-1);
                            
                            vgt_strain_rate = solver.vgt.strain(1/2 * [
                                0         gamma_dot;
                                gamma_dot 0]);
                            
                            vgt_stress_rate = C * vgt_strain_rate;
                            
                        case "Trusdell"
                            mat_L = solver.element.mat.v_half * B';
                            vgt_D = solver.vgt.B(B) * solver.element.vgt.v_half;
                            
                            vgt_stress_rate = C * vgt_D + solver.vgt.stress(mat_L * mat_stress + mat_stress * mat_L' - trace(mat_L) * mat_stress);
                            
                        otherwise
                            error("Constitutive model not implemented yet")
                    end
                    
                    k = 1;
                    vgt_stress(:, k) = vgt_stress + vgt_stress_rate * solver.dt(n-1);
                    
                    mat_strain_plastic(:, :, k) = zeros(2, 2);
                    strain_p_bar(k) = 0;
                    lambda(k) = 0;
                    mat_n_direction(:, :, k) = zeros(2, 2);
                    
                    phi(k) = plasticity.phi(deviatoric(vgt_stress(:, k)), strain_p_bar(k));
                    
                    % Plastic corrector
                    while(phi(k) > 0)
                        
                        % G = data.E / (2 * (1 + data.nu));
                        G = data.E / (1 - data.nu^2) * (1 - data.nu)/2;
                        
                        if(strain_p_bar(k) ~= 0)
                            H = 325e6 * 0.125 * (strain_p_bar(k))^(0.125-1);
                        else
                            H = 325e6;
                        end
                        
                        lambda(k) = phi(k) / (3*G + H);
                        
                        mat_n_direction(:, :, k) = plasticity.n_direction(deviatoric(vgt_stress(:, k)));
                        
                        mat_strain_plastic(:, :, k) = lambda(k) * mat_n_direction(:, :, k);
                        vgt_stress(:, k + 1) = vgt_stress(:, k) - C * solver.vgt.strain(mat_strain_plastic(:, :, k));
                        
                        strain_p_bar(k + 1) = plasticity.strain_p_bar(mat_strain_plastic(:, :, k));
                        
                        phi(k + 1) = plasticity.phi(deviatoric(vgt_stress(:, k + 1)), strain_p_bar(k + 1));
                        % phi(k + 1) = phi(k) + (vgt_n_direction(:, k)' * (- lambda(k) * C * vgt_n_direction(:, k)) - H * lambda(k));
                        
                        k = k + 1;
                    end
                    
                    solver.element.vgt.f_int = solver.element.vgt.f_int + weight_Xi * weight_Eta * solver.vgt.B(B)' * vgt_stress(:, k) * J;
                    solver.element.vgt.f_ext = solver.element.vgt.f_ext + weight_Xi * weight_Eta * solver.vgt.N(N)' * data.rho * [0; -9.81] * J; % Negligible with respect to the internal forces
                    
                    solver.element.stress(element_idx, i, j, :) = vgt_stress(:, k);
                    
                    if(i == 1 && j == 1 && element_idx == 1)
                        plotting_struct.sigma12(n) = vgt_stress(3, k);
                    end
                    
                end
            end
        end
        
        solver.element.f = solver.element.vgt.f_ext - solver.element.vgt.f_int;
        
        dt_critical = min(compute_min_distance(solver.element.mat.u) / sqrt(data.E / data.rho), dt_critical);
        
        solver.global.f(idx_nodes, n) = solver.global.f(idx_nodes, n) + solver.element.f;
        
    end
    
    solver.dt(n) = solver.dt(1);
    % solver.dt(n) = 0.89 * dt_critical;
    
    % 06) Compute acceleration
    solver.global.a(:, n) = solver.vgt.M \ solver.global.f(:, n);
    
    % 07) Nodal velocities at half-step
    if (n == 1)
        solver.global.v_half(:, n) = 0 + 1/2 * solver.global.a(:, n) * solver.dt(n);
    else
        solver.global.v_half(:, n) = solver.global.v_half(:, n-1) + solver.global.a(:, n) * solver.dt(n);
    end
    
    % 08) Enforce velocities BCs
    solver.global.v_half(boundary.velocities.index, n) = boundary.velocities.value;
    
    % 08) Check energy balance
    % solver.w_int(:, n) = solver.w_int(:, max(n-1, 1)) + (solver.u(:, n) - solver.u(:, max(n-1, 1)))' * ((solver.f(:, n) - solver.f(:, max(n-1, 1)))/2);
    % solver.w_ext(:, n) = solver.w_ext(:, max(n-1, 1)) + (solver.u(:, n) - solver.u(:, max(n-1, 1)))' * ((solver.f(:, n) - solver.f(:, max(n-1, 1)))/2);
    % solver.w_kin(:, n) = 1/2 * solver.v_half(:, n)' * solver.M * solver.v_half(:, n);
    % solver.w(n) = 1e-3 * (1e-3 * solver.w_kin(:, n) + solver.w_int(:, n) - solver.w_ext(:, n));
    % assert(abs(solver.w(:, n) - solver.w(:, max(n-1, 1))) < eps, "Energy Balance not verified")
    
    % 09) Update nodal displacements
    if (n == 1)
        solver.global.u(:, n) = solver.global.u(:, n) + solver.global.v_half(:, n) * solver.dt(n);
    else
        solver.global.u(:, n) = solver.global.u(:, n-1) + solver.global.v_half(:, n) * solver.dt(n);
    end
    
    % 10) Enforce displacements BCs
    solver.global.u(boundary.displacement.index, n) = boundary.displacement.value;
    
    % 11) Update time_step n and time t
    gamma(n+1) = atan(solver.global.u(3, n) / solver.global.u(4, n));
    plotting_struct.gamma(n) = gamma(n+1);
    
    t = t + solver.dt(n);
    n = n + 1;
    
end
toc

clear i j element_idx idx_nodes
clear n t C dt_critical
clear Xi Eta weight_Xi weight_Eta N B J vgt_stress mat_stress mat_L vgt_D vgt_stress_rate vgt_stress



%% Plots

if output.plots
    
    figure_results_analysis = figure('Name', 'Results analysis', 'NumberTitle', 'off');
    tiledlayout(1, 3);
    
    % Initial vs. Deformed configuration
    nexttile([1, 2]);
    hold on
    grid on
    axis equal;
    
    for element_idx = 1:size(mesh.connectivity, 1)
        
        connections = mesh.connectivity(element_idx, [1:end, 1]);
        idx_nodes = reshape((2 * connections.' + [-1, 0])', 1, []);
        
        u_mesh = mesh.coordinates(connections, :);
        u_deformed = solver.global.u(idx_nodes, end);
        
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
    
    plot(plotting_struct.gamma, plotting_struct.sigma12 * 1e-6, 'o', 'LineWidth', 2)
    
    title('Shear stress vs. Shear strain');
    xlabel('Shear strain [rad]');
    ylabel('Shear stress [MPa]');
    
end

clear element_idx connections idx_nodes u_mesh u_deformed


function S = deviatoric(vgt_stress)
stress = vgt_stress_to_mat(vgt_stress);
S = stress - trace(stress) / 3 * eye(size(stress, 1));
end
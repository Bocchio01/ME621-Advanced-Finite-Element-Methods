% ME621 - Advanced Finite Element Methods
% Assignment 2
% Tommaso Bocchietti
% A.Y. 2023/24 - W24

clc
clear variables
close all

%% Problem datas

rod = struct( ...
    'L0', 200, ...
    'A0', 300, ...
    'E', 70e+9, ...
    'rho0', 2700);

fem = struct( ...
    'N_elements', 5, ...
    'N_nodes_element', 4, ...
    'N_Gauss_points', 3, ...
    'stress_model', "FPK", ...
    'T_final', 0.01, ...
    'F_current', 0, ...
    'F_initial', 0, ...
    'F_final', 1000e3);

assert(fem.N_nodes_element >= 2, "The number of nodes per element must be at least 2")

N_nodes = (fem.N_nodes_element - 1) * fem.N_elements + 1;


%% Elemental matrices and vectors

el_handlers = struct( ...
    'N0', cell(1), ...
    'B0', cell(1), ...
    'f_int', NaN, ...
    'f_ext', NaN, ...
    'M', NaN);

shape_coefficients = compute_shape_coefficients_1D(fem.N_nodes_element);

N0_handlers = cell(1, fem.N_nodes_element);
B0_handlers = cell(1, fem.N_nodes_element);

for node_element_idx = 1:fem.N_nodes_element
    N0_handlers{node_element_idx} = @(Xi) polyval(shape_coefficients(node_element_idx, :), Xi);
    B0_handlers{node_element_idx} = @(Xi) polyval(polyder(shape_coefficients(node_element_idx, :)), Xi) * 2 / rod.L0;
end

el_handlers.N0 = @(Xi) cell2mat(cellfun(@(fun) fun(Xi), N0_handlers, 'UniformOutput', false));
el_handlers.B0 = @(Xi) cell2mat(cellfun(@(fun) fun(Xi), B0_handlers, 'UniformOutput', false));

el_handlers.f_int = @(rod, fem, u, B0) compute_f_int(rod, fem, u, B0);
el_handlers.f_ext = @(rod, fem, N0) compute_f_ext(rod, fem, N0);
el_handlers.M = compute_m(rod, fem, el_handlers.N0);

clear N0_handlers B0_handlers node_element_idx shape_coefficients


%% Algorithm setup

L = @(e) compute_gather_matrix(e, fem.N_nodes_element, N_nodes);

M_global = zeros(N_nodes);
for el_idx = 1:fem.N_elements
    M_global(:,:) = M_global(:,:) + L(el_idx)' * el_handlers.M * L(el_idx);
end

solver = struct( ...
    'M',  diag(sum(M_global, 2)), ...
    'dt', zeros(1, 1), ...
    'u', zeros(N_nodes, 1), ...
    'v_half', zeros(N_nodes, 1), ...
    'a', zeros(N_nodes, 1), ...
    'f', zeros(N_nodes, 1), ...
    'w', zeros(1, 1), ...
    'w_kin', zeros(1, 1), ...
    'w_int', zeros(1, 1), ...
    'w_ext', zeros(1, 1), ...
    'element', struct( ...
        'f', zeros(fem.N_nodes_element, 1), ...
        'f_int', zeros(fem.N_nodes_element, 1), ...
        'f_ext', zeros(fem.N_nodes_element, 1), ...
        'u', zeros(fem.N_nodes_element, 1), ...
        'v_half', zeros(fem.N_nodes_element, 1) ...
    ));

clear M_global el_idx


%% Main loop

n = 1;
t = 0;
eps = 1e-2;

while (t < fem.T_final)

    solver.f(:, n) = zeros(N_nodes, 1);
    dt_critical = +inf;

    % 05) Stress update algorithm
    for element_idx = 1:fem.N_elements

        solver.element.u = L(element_idx) * solver.u(:, max(n-1, 1));
        solver.element.v_half = L(element_idx) * solver.v_half(:, max(n-1, 1));

        solver.element.f_int = el_handlers.f_int(rod, fem, solver.element.u, el_handlers.B0);
        solver.element.f_ext = el_handlers.f_ext(rod, fem, el_handlers.N0);

        solver.element.f = solver.element.f_ext - solver.element.f_int;

        dt_critical = min(1e-3 * rod.L0 / (fem.N_elements) / (fem.N_nodes_element-1) / sqrt(rod.E / rod.rho0), dt_critical);

        solver.f(:, n) = solver.f(:, n) + L(element_idx)' * solver.element.f;

    end

    solver.dt(n) = 0.89 * dt_critical;

    % 06) Compute acceleration
    solver.a(:, n) = solver.M \ solver.f(:, n);

    % 07) Nodal velocities at half-step
    solver.v_half(:, n) = solver.v_half(:, max(n-1, 1)) + solver.a(:, n) * solver.dt(n);
    solver.v_half(1, n) = 0;

    % 09) Update nodal displacements
    solver.u(:, n) = solver.u(:, max(n-1, 1)) + solver.v_half(:, n) * solver.dt(n);

    % 10) Enforce displacements BCs
    solver.u(1, n) = 0;

    % 08) Check energy balance
    solver.w_int(:, n) = solver.w_int(:, max(n-1, 1)) + (solver.u(:, n) - solver.u(:, max(n-1, 1)))' * ((solver.f(:, n) - solver.f(:, max(n-1, 1)))/2);
    solver.w_ext(:, n) = solver.w_ext(:, max(n-1, 1)) + (solver.u(:, n) - solver.u(:, max(n-1, 1)))' * ((solver.f(:, n) - solver.f(:, max(n-1, 1)))/2);
    solver.w_kin(:, n) = 1/2 * solver.v_half(:, n)' * solver.M * solver.v_half(:, n);
    solver.w(n) = 1e-3 * (1e-3 * solver.w_kin(:, n) + solver.w_int(:, n) - solver.w_ext(:, n));
    assert(abs(solver.w(:, n) - solver.w(:, max(n-1, 1))) < eps, "Energy Balance not verified")

    % 11) Update time_step n and time t
    t = t + solver.dt(n);
    n = n + 1;

    fem.F_current = ((fem.F_final - fem.F_initial) / (fem.T_final - 0)) * t + fem.F_initial;

end

clear element_idx dt_critical



%% Plots

nexttile
hold on
grid on

plot(cumsum(solver.dt), 1e6 * solver.u(end, :), "LineWidth", 2);

title("End bar displacement");
xlabel("t [s]")
ylabel("u [\mum]")

legend(['F_{initial} =' num2str(fem.F_initial/1e3) ' kN'])

% % To observe the displacement of the nodes along the rod at various time
% % step. Given the high load at the end and the short period of time, we can
% % see that basically just the last two nodes starts moving, while the rest
% % of the bar remain in a steady condition.
% nexttile
% hold on
% grid on
% 
% for i = floor(linspace(1, size(solver.u, 2), 50))
%     plot(1:N_nodes, solver.u(:, i))
% end


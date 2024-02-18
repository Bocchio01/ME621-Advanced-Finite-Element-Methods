function m = compute_m(rod, fem, N0)

integrand_func = @(Xi) rod.rho0 * rod.A0 * N0(Xi)' * N0(Xi) * rod.L0 / 2;
m = method_gaussian_quadrature_integration(integrand_func, fem);

end
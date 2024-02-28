function f_ext = compute_f_ext(rod, fem, N0)

integrand_func = @(Xi) rod.rho0 * rod.A0 * 0 * N0(Xi)' * rod.L0 / 2;
f_body = method_gaussian_quadrature_integration(integrand_func, fem);

f_traction = N0(1)' * rod.A0 * fem.F_current - N0(-1)' * rod.A0 * fem.F_current;

f_ext = f_body + f_traction;

end
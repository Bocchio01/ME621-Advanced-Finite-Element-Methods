function [Kt, Fint] = model_taylor_1(U, geometry)

[L1, A1, E1, L2, A2, E2] = decompose_geometry(geometry);
[u, v] = decompose_u(U);

Kt = [[A1.*E1.*L1.^(-1),0];[0,A2.*E2.*L2.^(-1)]];

Fint = [A1.*E1.*L1.^(-1).*u;A2.*E2.*L2.^(-1).*v];

end

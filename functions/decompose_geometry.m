function [L1, A1, E1, L2, A2, E2] = decompose_geometry(geometry)

% fields = fieldnames(geometry);
%
% for i = 1:length(fields)
%     for j = 1:length(fields{i})
%         eval([fields{i}, num2str(j), ' = geometry.', fields{j}, '(',j,');']);
%     end
% end

L1 = geometry.L(1);
A1 = geometry.A(1);
E1 = geometry.E(1);

A2 = geometry.A(2);
L2 = geometry.L(2);
E2 = geometry.E(2);

end
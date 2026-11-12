% Load basis matrices as structure
matrices = load('pauli_basis_dim8.mat');

% Load each operator directly into memory
pauli_1x = matrices.pauli_1x;
pauli_1y = matrices.pauli_1y;
pauli_1z = matrices.pauli_1z;
pauli_2x = matrices.pauli_2x;
pauli_2y = matrices.pauli_2y;
pauli_2z = matrices.pauli_2z;
pauli_3x = matrices.pauli_3x;
pauli_3y = matrices.pauli_3y;
pauli_3z = matrices.pauli_3z;

pauli_1x_2x = matrices.pauli_1x_2x;
pauli_1x_2y = matrices.pauli_1x_2y;
pauli_1x_2z = matrices.pauli_1x_2z;
pauli_1y_2x = matrices.pauli_1y_2x;
pauli_1y_2y = matrices.pauli_1y_2y;
pauli_1y_2z = matrices.pauli_1y_2z;
pauli_1z_2x = matrices.pauli_1z_2x;
pauli_1z_2y = matrices.pauli_1z_2y;
pauli_1z_2z = matrices.pauli_1z_2z;

pauli_1x_3x = matrices.pauli_1x_3x;
pauli_1x_3y = matrices.pauli_1x_3y;
pauli_1x_3z = matrices.pauli_1x_3z;
pauli_1y_3x = matrices.pauli_1y_3x;
pauli_1y_3y = matrices.pauli_1y_3y;
pauli_1y_3z = matrices.pauli_1y_3z;
pauli_1z_3x = matrices.pauli_1z_3x;
pauli_1z_3y = matrices.pauli_1z_3y;
pauli_1z_3z = matrices.pauli_1z_3z;

pauli_2x_3x = matrices.pauli_2x_3x;
pauli_2x_3y = matrices.pauli_2x_3y;
pauli_2x_3z = matrices.pauli_2x_3z;
pauli_2y_3x = matrices.pauli_2y_3x;
pauli_2y_3y = matrices.pauli_2y_3y;
pauli_2y_3z = matrices.pauli_2y_3z;
pauli_2z_3x = matrices.pauli_2z_3x;
pauli_2z_3y = matrices.pauli_2z_3y;
pauli_2z_3z = matrices.pauli_2z_3z;

pauli_1x_2x_3x = matrices.pauli_1x_2x_3x;
pauli_1x_2x_3y = matrices.pauli_1x_2x_3y;
pauli_1x_2x_3z = matrices.pauli_1x_2x_3z;
pauli_1x_2y_3x = matrices.pauli_1x_2y_3x;
pauli_1x_2y_3y = matrices.pauli_1x_2y_3y;
pauli_1x_2y_3z = matrices.pauli_1x_2y_3z;
pauli_1x_2z_3x = matrices.pauli_1x_2z_3x;
pauli_1x_2z_3y = matrices.pauli_1x_2z_3y;
pauli_1x_2z_3z = matrices.pauli_1x_2z_3z;

pauli_1y_2x_3x = matrices.pauli_1y_2x_3x;
pauli_1y_2x_3y = matrices.pauli_1y_2x_3y;
pauli_1y_2x_3z = matrices.pauli_1y_2x_3z;
pauli_1y_2y_3x = matrices.pauli_1y_2y_3x;
pauli_1y_2y_3y = matrices.pauli_1y_2y_3y;
pauli_1y_2y_3z = matrices.pauli_1y_2y_3z;
pauli_1y_2z_3x = matrices.pauli_1y_2z_3x;
pauli_1y_2z_3y = matrices.pauli_1y_2z_3y;
pauli_1y_2z_3z = matrices.pauli_1y_2z_3z;

pauli_1z_2x_3x = matrices.pauli_1z_2x_3x;
pauli_1z_2x_3y = matrices.pauli_1z_2x_3y;
pauli_1z_2x_3z = matrices.pauli_1z_2x_3z;
pauli_1z_2y_3x = matrices.pauli_1z_2y_3x;
pauli_1z_2y_3y = matrices.pauli_1z_2y_3y;
pauli_1z_2y_3z = matrices.pauli_1z_2y_3z;
pauli_1z_2z_3x = matrices.pauli_1z_2z_3x;
pauli_1z_2z_3y = matrices.pauli_1z_2z_3y;
pauli_1z_2z_3z = matrices.pauli_1z_2z_3z;
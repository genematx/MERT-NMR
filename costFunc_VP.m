function result = costFunc_VP(yT, pars, c_ref, f_ref, t1, t2)
% Computes the value of the objective function for a 2D hypercomplex model
%  ----------------------------- Inputs: ---------------------------------
% yT - N_x_1 vector. Measured bicomplex data in the time domain.
% pars - K_x_4 matrix of model parameters. Each row corresponds to seprate
%        resonance, and different columns contain different parameters:
%        columns 1 and 2 - chemical shift in the first and second
%                          dimension;
%        columns 3 and 4 - peak relaxation rates (alpha) in the first and
%                          second dimension.
%        I.e., pars = [chsh, alpha].
% c_ref - 2_x_1 array of spectrometer frequencies in the first and second
%        dimensions in MHz
% f_ref - 2_x_1 array of frequency offsets in both dimensions in Hz
% t1, t2 - arrays of time samples in both dimesions (can have different
%        widths)
%  ----------------------------- Outputs: --------------------------------
% result - the value of the log-norm of the difference between the measured
%        data and the model specified by the array of parameters.

% Rearrange the parameters
pars = reshape(pars, [], 4);
chsh = pars(:, [1:2]);
alpha = pars(:, [3:4]);
yT = reshape(yT, [], 1);

% Compute the model matrix
Z = getModelMatrix(chsh, alpha, c_ref, f_ref, t1, t2);

% Find the value of the variable projection functional
ZZ = transpose(conj(Z))*Z;
Zy = transpose(conj(Z)) * yT;
iZZ = inv(ZZ);
result = real( transpose(conj(yT))*yT - transpose(conj(Zy)) * iZZ * Zy );

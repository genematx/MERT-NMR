function Z = getModelMatrix(chsh, alpha, c_ref, f_ref, t1, t2, varargin)
% Computes the bicomplex model matrix given the values of parameters.
%  ----------------------------- Inputs: ---------------------------------
% yT - N_x_1 vector. Measured bicomplex data in the time domain.
% chsh - K_x_2 matrix of model chemical shifts with rows corresponding to
%        separate resonances and each column containing chemical shifts in
%        both dimensions.
% alpha - K_x_2 matrix of relaxation rates that determine peak widths in
%        the spectrum; rows correspond to separate resonances and each
%        column contains parameters for direct and indirect dimensions.
% c_ref - 2_x_1 array of spectrometer frequencies in the first and second
%        dimensions in MHz
% f_ref - 2_x_1 array of frequency offsets in both dimensions in Hz
% t1, t2 - arrays of time samples in both dimesions (can have different
%        widths)
%
% Optional keyword parameters:
% 'tau' - 2_x_1 array of ringdow acquisition delays in microseconds
% 'sameAlpha' - boolean; default = False. If True, the peak widths of all
%        components in both dimensions will will be overwritten by the
%        widths of the first peak in the 'alpha' array.
% 'normalize' - boolean; default = True. Scale the model signals to make
%        them having norm=1.
%
%  ----------------------------- Output: ---------------------------------
% Z - N_x_K bicomplex matrix. The model matrix (defined in the time
%        domain). Each column of the matrix corresponds to a specific
%        resonance.

% Keep previously computed values
persistent chsh_old alpha_old Z_old tau_old

tau = [0 0];           % Default acquisition delays
normalize = true;
K = size(chsh, 1);     % Number of peaks
nt1 = numel(t1);
nt2 = numel(t2);
for i = 1 : 2 : numel(varargin)
    switch varargin{i}
        case 'tau'
            tau = varargin{i+1};
            tau = tau(:)';
        case 'sameAlpha'
            if varargin{i+1} == true
                alpha(:, 1) = alpha(1,1);
                alpha(:, 2) = alpha(1,2);
            end;
        case 'normalize'
            normalize = varargin{i+1};
    end;
end;

% Allocate memory for the model matrix
Z = bicomplex(NaN(nt1*nt2, K), NaN(nt1*nt2, K));

% Find which parameters have changed after the previous computation
if sum([numel(tau_old) == 0, numel(chsh) ~= numel(chsh_old), numel(alpha) ~= numel(alpha_old)])   % The first run
    indx_new = 1:K;        % Indices of the components that need to be recomputed
    Z = bicomplex(NaN(nt1*nt2, K), NaN(nt1*nt2, K));
elseif sum([tau_old ~= tau])
    indx_new = 1:K;        % Indices of the components that need to be recomputed
    Z = bicomplex(NaN(nt1*nt2, K), NaN(nt1*nt2, K));
else                       % Only (possibly some of) chemical shifts and/or alphas have changed
    [indx_new, ~] = find([chsh_old ~= chsh] + [alpha_old ~= alpha]);
    indx_new = unique(indx_new);
    Z = Z_old;
    tau = tau_old;
end;
indx_new = indx_new(:)';

% Compute Z for each peak separately
if numel(indx_new) > 0     % If something has changed
    omega = 2*pi*bsxfun(@minus, bsxfun(@times, chsh', c_ref), f_ref);    % An array of angular frequencies
    for i = indx_new
        S1 = exp(1i*omega(1, i)*(t1+tau(1)*(1e-06)) - alpha(i, 1)*t1) / sqrt(nt1);
        S2 = exp(1i*omega(2, i)*(t2'+tau(2)*(1e-06)) - alpha(i, 2)*t2') / sqrt(nt2);
        S1z = kron(real(S1), real(S2))+1i*kron(imag(S1), real(S2));
        S2z = kron(real(S1), imag(S2))+1i*kron(imag(S1), imag(S2));
        Z(:,i) = reshape(bicomplex(S1z, S2z), [], 1);
    end;
    
    % Normalize the model matrix
    if normalize
        Z = Z*diag(1./sqrt(diag(real(transpose(conj(Z))*Z))));
    end;
end;

% Keep the parameters
chsh_old = chsh;
alpha_old = alpha;
Z_old = Z;
tau_old = tau;


function Z = getModelMatrix(chsh, alpha, c_ref, f_ref, t1, t2, varargin)
% Computes the value of the objective function for a 2D hypercomplex model
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

% Compute Z for each peak separately
omega = 2*pi*bsxfun(@minus, bsxfun(@times, chsh', c_ref), f_ref);    % An array of angular frequencies
for i = 1 : K
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

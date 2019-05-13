function [result, Z, ZZ, Zy, a, PHS] = logPostHCph(yT, pars, c_ref, f_ref, t1, t2, varargin)
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
%
% Optional keyword parameters:
% 'tau' - 2_x_1 array of ringdow acquisition delays in microseconds
% 'mask' - an nf1_x_nf2 array 0 and 1. If suuplied, the comparison will be
%        carried out in the frequency domain. Only positions with mask==1
%        will be taken into account.
% 'noPhasing' - boolean. If True (default) no global phasing will be
%        applied to the data. Phases of each component will be estimated as
%        part of their bicomplex amplitudes.
% 'sameAlpha' - boolean; default = False. If True, the peak widths of all
%        components in both dimensions will will be overwritten by the
%        widths of the first peak in the 'pars' array.
% 'window' - window function to be applied in the time domain for
%        weighting. Choose from: 'none' (default), 'gauss', 'tukey',
%        'hann', 'rcos', 'qsine'.
%
%  ----------------------------- Outputs: --------------------------------
% result - the value of the log-norm of the difference between the measured
%        data and the model specified by the array of parameters.
% Z - N_x_K bicomplex matrix. The model matrix (defined in the time
%        domain). Each column of the matrix corresponds to a specific
%        resonance.

persistent chsh_old alpha_old Z_old ZZ_old Zy_old lPost_old tau_old;
noPhasing = 1;
pars = reshape(pars, [], 4);
chsh = pars(:, [1:2]);
alpha = pars(:, [3:4]);
tau = [0 0];           % Acquisition delays
K = size(chsh, 1);     % Number of peaks
nt1 = numel(t1);
nt2 = numel(t2);
sz = size(yT); n3 = sz(2);
mask = ones(1024, 1024);    % A mask for optimization in the frequency domain
wnd = 1;
for i = 1 : 2 : numel(varargin)
    switch varargin{i}
        case 'tau'
            tau = varargin{i+1};
            tau = tau(:)';
        case 'mask'    % ppm scales
            mask = varargin{i+1};
        case 'noPhasing'
            noPhasing = varargin{i+1};
        case 'sameAlpha'
            if varargin{i+1} == true;
                alpha(:, 1) = alpha(1,1);
                alpha(:, 2) = alpha(1,2);
            end;
        case 'window'
            % Define weighting windows in the time domains
            switch varargin{i+1}
                case 'none'
                    w1 = ones(nt1, 1);
                    w2 = ones(nt2, 1);
                case 'gauss'
                    w1 = window(@gausswin, nt1);
                    w2 = window(@gausswin, nt2);
                case 'tukey'
                    w1 = window(@tukeywin, nt1, 0.1);
                    w2 = window(@tukeywin, nt2, 0.1);
                case 'hann'
                    w1 = window(@hann, nt1);
                    w2 = window(@hann, nt2);     % Hann window
                case 'rcos'     % Raised cosine window
                    t = linspace(0, pi+pi/3, nt1)'; w1 = cos(t - pi/3)+1;
                    t = linspace(0, pi+pi/3, nt2)'; w2 = cos(t - pi/3)+1;
                case 'qsine'    % Raised and squared cosine window
                    t = linspace(0, pi+pi/3, nt1)'; w1 = ((cos(t - pi/3)+1)/2).^2;
                    t = linspace(0, pi+pi/3, nt2)'; w2 = ((cos(t - pi/3)+1)/2).^2;
                otherwise
                    error('Unknown window name.');
            end;
            wnd = w1(:) * w2(:)';
    end;
end;
if noPhasing; tau = [0 0]; end;

% Apply the window in time domain
yT = reshape(yT, nt1, nt2, n3);
for i = 1 : n3
    yT(:,:,i) = yT(:,:,i) .* wnd;
end;
yT = reshape(yT, [], n3);

% Find which parameters have changed from the previous computation
if sum([numel(tau_old) == 0, numel(chsh) ~= numel(chsh_old), numel(alpha) ~= numel(alpha_old)])   % The first run
    indx_new = 1:K;        % Indices of the components that need to be recomputed
    Z = bicomplex(NaN(nt1*nt2, K), NaN(nt1*nt2, K));
    ZZ = bicomplex(NaN(K, K), NaN(K, K));
    Zy = bicomplex(NaN(K, n3), NaN(K, n3));
elseif sum([tau_old ~= tau])
    indx_new = 1:K;        % Indices of the components that need to be recomputed
    Z = bicomplex(NaN(nt1*nt2, K), NaN(nt1*nt2, K));
    ZZ = bicomplex(NaN(K, K), NaN(K, K));
    Zy = bicomplex(NaN(K, n3), NaN(K, n3));
else                              % Only (possibly some of) chemical shifts and/or alphas have changed
    [indx_new, ~] = find([chsh_old ~= chsh] + [alpha_old ~= alpha]);
    indx_new = unique(indx_new);
    Z = Z_old;
    ZZ = ZZ_old;
    Zy = Zy_old;
    tau = tau_old;
    lPost = lPost_old;
end;
indx_new = indx_new(:)';

if numel(indx_new) > 0     % Something has changed
    
    % Form the inner product matrices and vectors
    omega = 2*pi*bsxfun(@minus, bsxfun(@times, chsh', c_ref), f_ref);
    for i = indx_new
        S1 = exp(1i*omega(1, i)*(t1+tau(1)*(1e-06)) - alpha(i, 1)*t1) / sqrt(nt1);
        S2 = exp(1i*omega(2, i)*(t2'+tau(2)*(1e-06)) - alpha(i, 2)*t2') / sqrt(nt2);
        S1z = kron(real(S1), real(S2))+1i*kron(imag(S1), real(S2));
        S2z = kron(real(S1), imag(S2))+1i*kron(imag(S1), imag(S2));
        Z(:,i) = reshape(bicomplex(S1z, S2z) .* wnd, [], 1);
    end;
    
    ZZ(indx_new, :) = transpose(conj(Z(:,indx_new))) * Z;
    ZZ(:, indx_new) = transpose(conj(ZZ(indx_new, :)));
    ZZ = (ZZ + transpose(conj(ZZ)))/2;
    Zy(indx_new, :) = transpose(conj(Z(:,indx_new))) * yT;
    iZZ = inv(ZZ); iZZ = (iZZ + transpose(conj(iZZ)))/2;
    yy = transpose(conj(yT))*yT;
    
    if noPhasing
        % Without phasing (assuming bicomplex amplitudes)
        a = iZZ * Zy;                         % a is bicomplex
        xThat = Z * a;
        PHS = 1;
        % lPost = -log(real(yy - transpose(conj(a))*ZZ*a));
    else
        % With phasing
        irZZ = inv(real(ZZ)); irZZ = (irZZ + transpose(conj(irZZ)))/2;
        % % Find the phasing term by solving the optimization problem
        min_opts = optimoptions('fminunc','Display', 'off', 'Algorithm','trust-region','GradObj','on');
        theta_hat = fminunc(@(x) fnc_phs(x, irZZ, Zy), [0, 0], min_opts);
        theta_hat = wrapToPi(theta_hat);
        PHS = bicomplex(exp(1i*theta_hat(1)), exp(1i*theta_hat(2)), 1);
        a = irZZ * real(Zy*conj(PHS));        % a is real
        if nnz(a < 0) > K/2;
            a = -a;
            PHS = -PHS;
        end;
        a = max(a, 0);
        xThat = Z * a * PHS;
    end;
    
    % Apply the mask and find the residual in the frequency domain
    if nnz(mask) ~= numel(mask)
        [nf1, nf2] = size(mask);
        indx_mask = find(mask);
        %     % yF = abs(modc(ffthc(reshape(yT, nt1, nt2), [nf1, nf2], 'hann', 'bicomplex')));
        %     % xFhat = abs(modc(ffthc(reshape(xThat, nt1, nt2), [nf1, nf2], 'hann', 'bicomplex')));
        yF = ffthc(reshape(yT, nt1, nt2), [nf1, nf2], 'rcos', 'bicomplex');
        xFhat = ffthc(reshape(xThat, nt1, nt2), [nf1, nf2], 'rcos', 'bicomplex');
        %     yF = ffthc(reshape(yT, nt1, nt2), [nf1, nf2], 'rcos', 'complex');
        %     xFhat = ffthc(reshape(xThat, nt1, nt2), [nf1, nf2], 'rcos', 'complex');
        r = (yF(indx_mask)-xFhat(indx_mask)) / sqrt(nnz(mask));
    else
        r = yT - xThat;
    end;
    
    % Find the norm of the residue
    lPost = -log(real(transpose(conj(r))*r));
    
end;

% Keep the parameters
chsh_old = chsh;
alpha_old = alpha;
Z_old = Z;
ZZ_old = ZZ;
Zy_old = Zy;
lPost_old = lPost;
tau_old = tau;
result = lPost;

function [f, g] = fnc_phs(x, irZZ, Zy)
    T = Zy*bicomplex(exp(-1i*x(1)), exp(-1i*x(2)), 1);
    f = - sum(diag(real(T)' * irZZ * real(T)));
    if nargout > 1
        g = [- sum(diag(imag1(T)' * irZZ * real(T))); ...
            - sum(diag(imag2(T)' * irZZ * real(T)))];
    end;


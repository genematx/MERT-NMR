clear all;
s = RandStream('mcg16807','Seed',0);
RandStream.setGlobalStream(s);

%% Default parameters
mode = 'bicomplex';   % 'complex';   % 
nf1 = 1024;        nf2 = 1024;       % Number of points in each dimension of the spectrum (after zero-filling)
indY = 1;

%% Choose an experiment to run
example = 1150;
switch example
    case 1150
        load(fullfile('_data', 'fitted_pars_1150_108_ph'));
        true_chsh = best_chsh; true_alpha = best_alpha;
        fid = fopen(fullfile('data', '1150my_t1t2.spc'));
        t_echo = [0.08 0.16 0.24 0.4 0.56 0.64 0.72 0.8];
    case 1151
        load(fullfile('_data', 'fitted_pars_1150_108_ph')); true_chsh = best_chsh; true_alpha = best_alpha;
        fid = fopen(fullfile('data', '1151my_t1t2.spc'));
        t_echo = [0.001 0.005 0.01 0.018 0.035 0.055 0.08 0.11 0.15 0.2];
    case 1155
        load(fullfile('_data', 'fitted_pars_1150_108_ph')); true_chsh = best_chsh; true_alpha = best_alpha;
        fid = fopen(fullfile('data', '1155my_t1t2.spc'));
        t_echo = [0.0144 0.0288 0.0432 0.0576 0.072 0.0864 0.1008 0.1152 0.1296 0.144];
end;
best_chsh = true_chsh; best_alpha = true_alpha;% / pi

%% Load the data, define parameters
% Experiments 115x
if and(example > 1100, example < 1200)
%% Load the data
procdat = fread(fid,inf,'float',0,'l');
fclose(fid);
dim1 = 1024; %direct dimension points (we get these for "free")
dim2 = 256;  %indirect dimension points (these we would want to undersample)
dim3 = numel(t_echo);    %time delays for relaxation fit
block1 = 128;
block2 = 16;
block3 = 2;

nb1 = ceil(dim1/block1);%floor(length(procdat)/dim2/block1);
nb2 = ceil(dim2/block2);%floor(length(procdat)/(nb1*block1))/block2);
nb3 = ceil(dim3/block3);%floor(length(procdat)/(nb1*block1))/block2);
temp = reshape(procdat,[block1,block2,block3,nb1,nb2,nb3]);
temp = permute(temp,[1 4 2 5 3 6]);
temp = reshape(temp,[nb1*block1,nb2*block2,nb3*block3]);

yT = bicomplex((temp(1:2:end, 2:2:end, :)+1i*temp(2:2:end, 2:2:end, :)), ...
    -(temp(1:2:end, 1:2:end, :)+1i*temp(2:2:end, 1:2:end, :)));
yT = yT(1:440, :, :);               % Cut the part left due to the initial delay
yT = yT / norm(norm(yT(:,:,2)));    % Just scale everything by a constant
%% Define the signal parameters
nt = size(yT); nt1 = nt(1); nt2 = nt(2);      % Number of points in each time dimension
% nf1 = nt1; nf2 = nt2;       % Number of points in each dimension of the spectrum
c1_ref = 500.13;                    c2_ref = 50.677748;
f1_ref = 2360;                      f2_ref = 6081.32976;
f1_ref = f1_ref + 0.038*c1_ref;     f2_ref = f2_ref - 2.8*c2_ref;     % Correction to make the scales match with the provided analyzed data
dt1 = 1e-04;                        dt2 = 6.32e-04;
fs1 = 1/dt1;                        fs2 = 1/dt2;
f1 = linspace(-fs1/2, fs1/2, nf1+1); f1 = f1(1:end-1);
f2 = linspace(-fs2/2, fs2/2, nf2+1); f2 = f2(1:end-1);
del1 = (f1 + f1_ref ) /c1_ref ;  %   Frequency scale in ppm (for even number of points)
del2 = (f2 + f2_ref ) /c2_ref ;  % 
t1 = [0:nt1-1]'*dt1;        % Time in seconds
t2 = [0:nt2-1]'*dt2;
n_echo = numel(t_echo);
%% Cut only an interesting part in the spectrum (i.e. 6.2...10.2 ppm in the 1H direction)
ind_cut = [234:361];      % Indices in the frequency domain to cut (an even number)          ind_cut = [80:207];
yt1 = NaN(nt1, 2*nt2, n_echo);     % Data in time domain, complex values along the 1st dimension, R-I-R-I-... along the second
yt1(:, 1:2:end, :) = cmpl1(yT(:,:,:)); yt1(:, 2:2:end, :) = cmpl2(yT(:,:,:));
yf1t2 = fftshift(fft(yt1, nt1, 1), 1);
yf1t2 = yf1t2(ind_cut, :, :);
f1 = linspace(-fs1/2, fs1/2, nt1+1); f1 = f1(1:end-1);  % Frequency vector for the computed FT
yt1t2 = ifft(ifftshift(yf1t2, 1), [], 1);   % Back to the time domain
yT = bicomplex(yt1t2(:, 1:2:end, :), yt1t2(:, 2:2:end, :));

% Redefine the parameters for the first dimension
nt1 = numel(ind_cut);      % New number of time points
c1_ref = c1_ref;   % Stays the same
f1 = f1(ind_cut) + f1_ref;
f1_ref = f1(nt1/2+1);   % Set the middle point in the spectrum to be the reference frequency
f1 = f1 - f1_ref;  % New frequency range (in Hz)
fs1 = 2*abs(f1(1));     % Sampling frequency
dt1 = 1/fs1;       % Dwell time
% Generate new range vector
f1 = linspace(-fs1/2, fs1/2, nf1+1); f1 = f1(1:end-1);
del1 = (f1 + f1_ref) /c1_ref ;  %   Frequency scale in ppm (for even number of points)

t1 = [0:nt1-1]'*dt1;        % Time in seconds
[tx, ty] = meshgrid(t2, t1);
[delx, dely] = meshgrid(del1, del2);
[fx, fy] = meshgrid(f1, f2);   % Frequency grid
Dgrd = [delx(:) dely(:)];
end;

%% Main optimization
n_peaks = size(true_chsh, 1);    % Total number of peaks
c_ref = [c1_ref; c2_ref];
f_ref = [f1_ref; f2_ref];
omega = 2*pi*bsxfun(@minus, bsxfun(@times, true_chsh', c_ref), f_ref);
alpha = true_alpha';
if strcmp(mode, 'bicomplex')
%% Construct the model signals - HC
[result, Z, ZZ, Zy] = logPostHC(reshape(yT(:, :, 1), [], 1), [true_chsh, true_alpha], c_ref, f_ref, t1, t2);
ZZi = inv(ZZ);
b = ZZi*Zy;     % Inferred b
Zb = Z * b;
xT = reshape(Zb, nt1, nt2);                     % Modeled signal
% xF = fftshift(fftshift(fft2(xT, nf1, nf2), 1), 2) / sqrt(nt1*nt2);      % frequency-domain representation of x
%% Rescale the data to fit the model (linear scaling and offset)
Y = [reshape(yT(:,:,1), [], 1) bicomplex(ones(nt1*nt2,1), zeros(nt1*nt2, 1))];    % Observed data
YY = transpose(conj(Y)) * Y;
iYY = inv(YY);
p = iYY * (transpose(conj(Y))*reshape(xT, [], 1));        % Scaling coefficient and an offset
yT = yT .* p(1) + p(2);
clear p Y YY iYY;
% return;
%% Fit the model closer to the data
yT1 = reshape(yT(:,:,1), [], 1);    % The first panel; will be used to fit the signal parameters
yyT1 = transpose(conj(yT1))*yT1;

% Initialization of the parameters
best_pars = [best_chsh best_alpha];
% best_pars(:, 3:4) = 1;
allPtcl = best_pars(:);
% Set the range limits for optimization
D = pdist2(true_chsh, true_chsh); Dmin = D; Dmin(1:n_peaks+1:end) = Inf;
Dmin = min(Dmin)'; Dmin = min(Dmin, 0.5);   % Distance to the closest peak
lwr_pars = bsxfun(@minus, best_pars, [min(Dmin, 0.05)/2 min(Dmin, 0.5)/2 35*ones(n_peaks,1) 30*ones(n_peaks,1)]);
lwr_pars(:, 3) = 0; lwr_pars(:, 4) = 0;
upr_pars = bsxfun(@plus, best_pars, [min(Dmin, 0.05)/2 min(Dmin, 0.5)/2 20*ones(n_peaks,1) 20*ones(n_peaks,1)]);
upr_pars(:, 3) = 200; upr_pars(:, 4) = 200;
% Evaluate the log Posterior
[allPost, Z, ZZ, Zy] = logPostHC(yT1, allPtcl(:, end)', c_ref, f_ref, t1, t2);
clear logPostHC;

% % Fit the frequencies iteratively
% indx_optm = reshape(1:(n_peaks*4), [], 4);
% indx_optm = indx_optm(:, 3:4);
% indx_optm = indx_optm';
% indx_optm = indx_optm(:)';
% indx_optm = mat2cell(indx_optm, 1, 2*ones(1, n_peaks));
% indx_optm = indx_optm([randperm(numel(indx_optm)) randperm(numel(indx_optm)) randperm(numel(indx_optm))]);
% % % Fit peak by peak
% % indx_optm = reshape(1:(n_peaks*4), [], 4);
% % indx_optm = indx_optm([203], 1:4);
% % indx_optm = indx_optm';
% % indx_optm = indx_optm(:)';
% % indx_optm = mat2cell(indx_optm, 1, 4*ones(1, numel(indx_optm)/4));

% % % Fit tau's (values are in microseconds)
allTau = best_tau(:);
logPost_tau = @(x) -logPostHCph(yT1, allPtcl(:, end), c_ref, f_ref, t1, t2, 'tau', x, 'noPhasing', 0);
[new_tau, mcPost, exitflag, output] = fmincon(logPost_tau, ...
    allTau(:, end)',[],[],[],[],[-100, -100], [100, 100], [], optimset('Display', 'iter'));
allPtcl(:, end+1) = allPtcl(:, end);    % Keep all results
allPost(end+1) = mcPost;
allTau(:, end+1) = new_tau(:);

sameAlpha = false;
numTryPts = 0;
% With built-in phasing
for kk = 101   % [4     9    10    18    20    39    45    48    49    71    73   100   119   121   123   152   153   154   172   175   207] % randperm(n_peaks)    %  numel(indx_optm)
    if ~mod(kk, 5); clear logPostHCph logPost_rdcd; end;
    indx_peak = kk;
    Dngb = pdist2(best_chsh, best_chsh(indx_peak, :), 'mahalanobis');    % Find neighbours
    indx_nghb = find(Dngb < 0.13);
    indx_optm{kk} = bsxfun(@plus, [0,1,2, 3]*n_peaks, indx_nghb);
    indx_optm{kk} = indx_optm{kk}(:)';
    
    pars_new = allPtcl(:, end)';
    fprintf('Analyzing step #%d of %d.\n', kk, numel(indx_optm));
    if sameAlpha
        pars_new((2*n_peaks+1):(3*n_peaks)) = pars_new(2*n_peaks+1);
        pars_new((3*n_peaks+1):(4*n_peaks)) = pars_new(3*n_peaks+1);
    end;

% ----- Create a mask for optimization in frequency domain
    [indx_peak, dummy] = ind2sub([n_peaks, 4], indx_optm{kk});
    indx_peak = unique(indx_peak);       % Only peaks that are being optimized now
    distGrd = pdist2(Dgrd, best_chsh(indx_peak, :), 'mahalanobis');
    distGrd = min(distGrd, [], 2);
    mask = zeros(nf1, nf2);
    mask(distGrd < 0.16) = 1;
    mask = fliplr(rot90(mask, -1));
%     mask = 1;
    
% ----- Form parameter-selecting matrices
    K = ones(1, 4*n_peaks); K(indx_optm{kk}) = 0;
    K = diag(K);   % Keep
    if nnz(K) == 0; K = 0; O = 1;
    else
        O = eye(4*n_peaks) - K; O = O(indx_optm{kk}, :);   % Optimize
        K = pars_new*K;
    end;
    logPost_rdcd = @(x) -logPostHCph(yT1, x*O+K, c_ref, f_ref, t1, t2, 'tau', allTau(:, end), 'mask', mask, 'noPhasing', 1);
    pars_0 = [pars_new(indx_optm{kk})];
   
% --- MATLAB global search
    opts = optimoptions(@fmincon,'Algorithm','sqp');
    if numTryPts == 0
        [gsPtcl, gsPost, flg, output] = fmincon(@(x) logPost_rdcd(x), pars_0, ...
            [],[],[],[],lwr_pars(indx_optm{kk}),upr_pars(indx_optm{kk}), [], opts);
    else
        gsProblem = createOptimProblem('fmincon','objective',@(x) logPost_rdcd(x),...
            'x0', pars_0, 'lb', [lwr_pars(indx_optm{kk})], 'ub', [upr_pars(indx_optm{kk})], 'options', opts);
        gs = GlobalSearch('Display', 'iter', 'NumTrialPoints', 10, 'NumStageOnePoints', 5, 'DistanceThresholdFactor', 0.9);    % gs = GlobalSearch('Display', 'iter', 'NumTrialPoints', 250, 'NumStageOnePoints', 150, 'DistanceThresholdFactor', 0.9);
        [gsPtcl, gsPost, flg, og, solns] = run(gs,gsProblem);
    end;
        
    pars_0_all = gsPtcl*O + K;
    allPtcl(:, end+1) = pars_0_all';    % Keep all results
    allPost(end+1) = gsPost;
    
    best_alpha = allPtcl(2*n_peaks+1:4*n_peaks, end);
    best_alpha = reshape(best_alpha, [], 2);
    best_chsh = allPtcl(1:2*n_peaks, end);
    best_chsh = reshape(best_chsh, [], 2);
    save(['fitting_progress_' num2str(example) '_bicmp'], 'allTau', 'allPtcl', 'allPost', 'best_alpha', 'best_chsh', 'upr_pars', 'lwr_pars');    % 
%     if and(~best_new, ii >= 10); break; end;
    fprintf('\n');
end;
best_pars = reshape(allPtcl(:, end), [], 4);

%% Compute spectra for plotting
[result, Z, ZZ, Zy, a, PHS] = logPostHCph(reshape(yT(:, :, 1), [], 1), allPtcl(:, end), c_ref, f_ref, t1, t2, 'tau', allTau(:, end), 'noPhasing', 1);
xT = Z * a * PHS;
xT = reshape(xT, nt1, nt2);                     % Modeled signal

% % Absolute values of bicomplex data
% peaks_plot_1 = abs(modc(ffthc(yT(:,:,indY), [nf1, nf2], 'none', 'bicomplex')));
% peaks_plot_2 = abs(modc(ffthc(xT(:,:), [nf1, nf2], 'none', 'bicomplex')));

% Real spectra of reduced complex data
peaks_plot_1 = real(ffthc(yT(:,:,indY), [nf1, nf2], 'rcos', 'complex'));
peaks_plot_2 = real(ffthc(xT(:,:), [nf1, nf2], 'rcos', 'complex'));

elseif strcmp(mode, 'complex')
%% Convert the bicomplex signals to the complex form
yt1 = NaN(nt1, 2*nt2, n_echo);     % Data in time domain, complex values along the 1st dimension, R-I-R-I-... along the second
yt1(:, 1:2:end, :) = cmpl1(yT); yt1(:, 2:2:end, :) = cmpl2(yT);
yf1t2 = fftshift(fft(yt1, nt1, 1), 1);
yf1t2 = real(yf1t2(:, 1:2:end, :)) + 1i*real(yf1t2(:, 2:2:end, :));
yTc = ifft(ifftshift(yf1t2, 1), nt1, 1);     % yTc = fftshift(yTc, 1);        % Move the sample corresponding to t=0 to the center of array
t1 = ([-floor(nt1/2):floor(nt1/2)])'*dt1; t1 = t1(1:end-1);       % Time in seconds
t1 = ifftshift(t1);

[tx, ty] = meshgrid(t2, t1);
[delx, dely] = meshgrid(del2, del1);

% %% Convert the bicomplex signals to the complex form -- 2
% for i = 1 : n_echo
%     temp = ffthc(yT(:,:,i), [nt1, nt2], 'rcos', 'bicomplex');
%     yFc(:,:,i) = real(temp) + 1i*imag12(temp);
%     yTc(:,:,i) = ifft2(ifftshift(ifftshift(yFc(:,:,i), 1), 2));
% end;
% t1 = ([-floor(nt1/2):floor(nt1/2)])'*dt1; t1 = t1(1:end-1);       % Time in seconds
% t2 = ([-floor(nt2/2):floor(nt2/2)])'*dt2; t2 = t2(1:end-1);       % Time in seconds
% t1 = ifftshift(t1);
% t2 = ifftshift(t2);
% 
% [tx, ty] = meshgrid(t2, t1);
% [fx, fy] = meshgrid(del2, del1);
%% Construct the model signals -- Complex
Z = NaN(nt1, nt2, n_peaks);
ZZ = NaN(n_peaks, n_peaks);
for i = 1 : n_peaks
    Z(:,:,i) = kron(exp(1i*omega(1, i)*t1 - true_alpha(i, 1)*abs(t1)), exp(1i*omega(2, i)*t2' - true_alpha(i, 2)*abs(t2')));
    Z(:,:,i) = Z(:,:,i)/norm(Z(:,:,i));
    for j = 1 : i
        ZZ(i, j) = reshape(Z(:,:,i), [], 1)' * reshape(Z(:,:,j), [], 1);
        ZZ(j, i) = ZZ(i, j)';
    end;
end;
ZZ = (ZZ + ZZ')/2;   % Make sure ZZ is Hermitian to avoid numerical errors
% intn = true_intn'; b = transpose( intn .* sum(exp(1i*(omega(:)'*tau + theta(:)'))) );    % Modeled b
Zy = reshape(Z, [], n_peaks)' * reshape(yTc(:, :, 2), [], 1);
b = inv(ZZ)*Zy; b = b / norm(b);     % Inferred b
% b = ones(K, 1);
Zb = bsxfun(@times, Z, reshape(b, 1, 1, []));
xTc = sum(Zb, 3);
%% Rescale the data to fit the model (linear scaling and offset)
Y = [reshape(yTc(:,:,1), [], 1) ones(nt1*nt2,1)];    % Observed data
YY = transpose(conj(Y)) * Y;
iYY = inv(YY);
p = iYY * (transpose(conj(Y))*reshape(xTc, [], 1));        % Scaling coefficient and an offset
yTc = yTc .* p(1) + p(2);
clear p Y YY iYY;
% return;
%% Fit the model closer to the data
yT1 = yTc(:,:,1);    % The first panel; will be used to fit the signal parameters
% Defina a prior for the amplitudes (Gaussian)
S0 = eye(n_peaks);
m0 = zeros(n_peaks, 1);
delta = 1e+42;
iS0 = inv(S0)/(delta^2);
% Define a prior for the variance of noise (Inverse Gamma)
alph_sigma2_prior = 2;
beta_sigma2_prior = 0.1;         % Prior on sigma^2~IG(alph_sigma2_prior, beta_sigma2_prior)

% Initialization of the parameters
best_pars = [true_chsh true_alpha];
allPtcl = best_pars(:);
% Set the range limits for optimization
D = pdist2(true_chsh, true_chsh); Dmin = D; Dmin(1:n_peaks+1:end) = Inf;
Dmin = min(Dmin)'; Dmin = min(Dmin, 0.5);   % Distance to the closest peak
lwr_pars = bsxfun(@minus, best_pars, [Dmin/3 Dmin/3 35*ones(n_peaks,1) 30*ones(n_peaks,1)]);
lwr_pars(:, 3) = 1; lwr_pars(:, 4) = 1;
upr_pars = bsxfun(@plus, best_pars, [Dmin/3 Dmin/3 20*ones(n_peaks,1) 20*ones(n_peaks,1)]);
upr_pars(:, 3) = 150; upr_pars(:, 4) = 150;
% Evaluate the log Posterior
omega = 2*pi*bsxfun(@minus, bsxfun(@times, best_pars(:, 1:2)', c_ref), f_ref);
alpha = best_pars(:, 3:4)';
Z = NaN(nt1, nt2, n_peaks);
ZZ = NaN(n_peaks, n_peaks);
for i = 1 : n_peaks
    Z(:,:,i) = kron(exp(1i*omega(1, i)*t1 - alpha(1, i)*abs(t1)), exp(1i*omega(2, i)*t2' - alpha(2, i)*abs(t2')));
    Z(:,:,i) = Z(:,:,i)/norm(Z(:,:,i));
    for j = 1 : i
        ZZ(i, j) = reshape(Z(:,:,i), [], 1)' * reshape(Z(:,:,j), [], 1);
        ZZ(j, i) = ZZ(i, j)';
    end;
end;
ZZ = (ZZ + ZZ')/2;   % Make sure ZZ is Hermitian to avoid numerical errors
Zy = reshape(Z, [], n_peaks)' * yT1(:);
iSc = iS0 + ZZ;
if rank(iSc) < n_peaks
    warning('The rank of ZZ is less than the number of peaks.\n');
    iSc = iSc + 0.00000001*eye(n_peaks);
end;
Sc = inv(iSc); Sc = (Sc + Sc')/2;
mc = Sc*(Zy+iS0*m0);
% Compute the log posterior
lPost = log(det(Sc)) - (alph_sigma2_prior + nt1*nt2)*log(beta_sigma2_prior + real(yT1(:)'*yT1(:) + m0'*iS0*m0 - mc'*iSc*mc));        % Terms that depend on many variables
lPost = lPost - nt1*nt2*log(pi);    % Terms that depend only on sigma and/or delta
lPost = real(lPost);    % Make sure it is real

allPost = lPost;
best_Z = Z;
best_ZZ = ZZ;
best_Zy = Zy;

% Set distributions to sample the parameters
for j = 1 : numel(best_pars)
        mu_gmm = (upr_pars(j) + lwr_pars(j)) / 2;
        sigma_gmm = (upr_pars(j) - lwr_pars(j)) / 5;
        p_gmm = 1;
        gmm{1, j} = gmdistribution(mu_gmm, sigma_gmm, p_gmm);
end;

% Update parameters iteratively
maxiter_IS = 100;
n_ptcl = 150;       % Number of random particles to sample
p_unif = 0.5;       % The ratio of particles sampled uniformly
indx_optim = reshape([1:4*n_peaks], n_peaks, 4);
% indx_optim = indx_optim(problematic_peaks, 3:4);   % Choose which parameters to optimize
indx_optim = indx_optim(:, 1:4);   % Choose which parameters to optimize
indx_optim = indx_optim(:)';
for ii = 2 : maxiter_IS
    best_new = false;
    fprintf('Iteration # %d. Parameter #', ii);
    for j = indx_optim(randperm(numel(indx_optim)))    % randperm(4*K)
        pars_new = allPtcl(:, end);
        i = mod(j-1, n_peaks) + 1;     % Index of the peak affected by the j'th parameter; will change only the corresponding rows/columns of Z nad ZZ
        if j > 2*n_peaks; p_unif = 0.5; else p_unif = 0.1; end;
        fprintf('%d ', j);
        % Generate new samples
        particles = repmat(pars_new, 1, n_ptcl);
        particles(j, :) = NaN;
        % Sample some particles uniformly
        particles(j, 1:round(p_unif*n_ptcl)) = bsxfun(@plus, bsxfun(@times, rand(round(p_unif*n_ptcl), 1), upr_pars(j)-lwr_pars(j)), lwr_pars(j));
        % Sample the rest of new particles and make sure they are in the allowed range
        indx_nan = find(isnan(particles(j, :)));
        while numel(indx_nan) > 0
            ptcl_new = random(gmm{ii-1, j}, numel(indx_nan));
            ptcl_new(or(ptcl_new<lwr_pars(j), ptcl_new>upr_pars(j))) = NaN;
            particles(j, indx_nan) = ptcl_new;
            indx_nan = find(isnan(particles(j, :)));
        end;
        W = -log((1-p_unif)*pdf(gmm{ii-1, j}, particles(j, :)') + p_unif);    % Initial distribution from which the new particles were sampled
                
        % Compute the posterior
        Post_ptcl = NaN(n_ptcl, 1);
        Z = best_Z; ZZ = best_ZZ; Zy = best_Zy;
        best_post_new = allPost(end);
        for p = 1 : n_ptcl
            crnt_ptcl = reshape(particles(:, p), n_peaks, 4);
            omega = 2*pi*bsxfun(@minus, bsxfun(@times, crnt_ptcl(i, 1:2)', c_ref), f_ref);
            alpha = crnt_ptcl(i, 3:4)';
            Z(:,:,i) = kron(exp(1i*omega(1)*t1 - alpha(1)*abs(t1)), exp(1i*omega(2)*t2' - alpha(2)*abs(t2')));
            Z(:,:,i) = Z(:,:,i)/norm(Z(:,:,i));
            ZZ(i, :) = reshape(Z(:,:,i), [], 1)' * reshape(Z, [], n_peaks);
            ZZ(:, i) = ZZ(i, :)'; 
            ZZ = (ZZ + ZZ')/2;   % Make sure ZZ is Hermitian to avoid numerical errors
            Zy(i) = reshape(Z(:,:,i), [], 1)' * yT1(:);
            iSc = iS0 + ZZ;
            if rank(iSc) < n_peaks
                warning('The rank of ZZ is less than the number of peaks.\n');
                iSc = iSc + 0.00000001*eye(n_peaks);
            end;
            Sc = inv(iSc); Sc = (Sc + Sc')/2;
            mc = Sc*(Zy+iS0*m0);
            % Compute the log posterior
            lPost = log(det(Sc)) - (alph_sigma2_prior + nt1*nt2)*log(beta_sigma2_prior + real(yT1(:)'*yT1(:) + m0'*iS0*m0 - mc'*iSc*mc));        % Terms that depend on many variables
            lPost = lPost - nt1*nt2*log(pi);    % Terms that depend only on sigma and/or delta
            lPost = real(lPost);    % Make sure it is real
            Post_ptcl(p) = lPost;
            
            % Keep the best matrices Z
            if lPost > best_post_new
                best_post_new = lPost;
                best_Z_new = Z;
                best_ZZ_new = ZZ;
                best_Zy_new = Zy;
                best_ptcl_new = crnt_ptcl(:);
            end;
        end;
        
        % Update the weights
        W = W + Post_ptcl(:); 
        W = bsxfun(@minus, W, max(W));
        W = exp(W);
        W = bsxfun(@rdivide, W, sum(W));     % Make all weights sum to 1 (define a categorical distribution)
        
        % Generate a new GMM
        [vals indx_gmm] = sort(W, 'descend');
        vals = cumsum(vals);
        indx_gmm = indx_gmm(1:max(1, min(find(vals > 0.97, 1, 'first'), 3) ));    % Do not take too many samples (3 is enough)
        mu_gmm = particles(j, indx_gmm)';
        sigma_gmm = (3/4*n_ptcl)^(-1/5) * std(particles(j, :));
        p_gmm = W(indx_gmm);
        gmm{ii, j} = gmdistribution(mu_gmm, sigma_gmm, p_gmm);
        
        % Find the MAP estimate of the new parameter
        xx = linspace(lwr_pars(j), upr_pars(j), 10000)';
        yy = pdf(gmm{ii, j}, xx);
        [val ind] = max(yy);
        
        % Try the MAP estimate
        pars_new(j) = xx(ind);
        crnt_ptcl = reshape(pars_new, n_peaks, 4);
        omega = 2*pi*bsxfun(@minus, bsxfun(@times, crnt_ptcl(i, 1:2)', c_ref), f_ref);
        alpha = crnt_ptcl(i, 3:4)';
        Z(:,:,i) = kron(exp(1i*omega(1)*t1 - alpha(1)*abs(t1)), exp(1i*omega(2)*t2' - alpha(2)*abs(t2')));
        Z(:,:,i) = Z(:,:,i)/norm(Z(:,:,i));
        ZZ(i, :) = reshape(Z(:,:,i), [], 1)' * reshape(Z, [], n_peaks);
        ZZ(:, i) = ZZ(i, :)';
        ZZ = (ZZ + ZZ')/2;   % Make sure ZZ is Hermitian to avoid numerical errors
        Zy(i) = reshape(Z(:,:,i), [], 1)' * yT1(:);
        iSc = iS0 + ZZ;
        if rank(iSc) < n_peaks
            warning('The rank of ZZ is less than the number of peaks.\n');
            iSc = iSc + 0.00000001*eye(n_peaks);
        end;
        Sc = inv(iSc); Sc = (Sc + Sc')/2;
        mc = Sc*(Zy+iS0*m0);
        % Compute the log posterior
        lPost = log(det(Sc)) - (alph_sigma2_prior + nt1*nt2)*log(beta_sigma2_prior + real(yT1(:)'*yT1(:) + m0'*iS0*m0 - mc'*iSc*mc));        % Terms that depend on many variables
        lPost = lPost - nt1*nt2*log(pi);    % Terms that depend only on sigma and/or delta
        lPost = real(lPost);    % Make sure it is real
        
        % Keep the best matrices Z
        if lPost > best_post_new
            best_post_new = lPost;
            best_Z_new = Z;
            best_ZZ_new = ZZ;
            best_Zy_new = Zy;
            best_ptcl_new = crnt_ptcl(:);
        end;
        
        % Update the best particle
        if best_post_new > allPost(end)
            allPost(end+1) = best_post_new;
            allPtcl(:, end+1) = best_ptcl_new;
            best_Z = best_Z_new;
            best_ZZ = best_ZZ_new;
            best_Zy = best_Zy_new;
            best_new = true;
        end;
        
    end;
    
    
    best_alpha = allPtcl(2*n_peaks+1:4*n_peaks, end);
    best_alpha = reshape(best_alpha, [], 2);
    best_chsh = allPtcl(1:2*n_peaks, end);
    best_chsh = reshape(best_chsh, [], 2);
    save(['fitting_progress_' num2str(example) '_cmplx'], 'best_pars_all', 'best_post', 'best_alpha', 'best_chsh');
    if and(~best_new, ii >= 10); break; end;
    fprintf('\n');
end;
    %% Compute spectra for plotting
    nt1half = ceil(numel(t1)/2);    % Number of the first-dimension samples in the first half
    yT_zf = [yTc(1:nt1half, :, :); zeros(nf1-nt1, nt2, n_echo); yTc(nt1half+1:end, :, :)];   % Zero-filling
    xT_zf = [xTc(1:nt1half, :); zeros(nf1-nt1, nt2); xTc(nt1half+1:end, :)];   % Zero-filling
    % Compute the FT
    yF = fftshift(fftshift(fft2(yT_zf, nf1, nf2), 1), 2) / sqrt(nt1*nt2);      % frequency-domain representation of y
    xF = fftshift(fftshift(fft2(xT_zf, nf1, nf2), 1), 2) / sqrt(nt1*nt2);      % frequency-domain representation of y
    peaks_plot_1 = real(yF(:,:,indY));
    peaks_plot_2 = real(xF);
end;

%% Compute spectra for plotting; if needed convert from bicomplex to complex form
% if strcmp(mode, 'complex')
%     %% Convert the bicomplex signals to the complex form
%     yt1 = NaN(nt1, 2*nt2, n_echo);     % Data in time domain, complex values along the 1st dimension, R-I-R-I-... along the second
%     yt1(:, 1:2:end, :) = cmpl1(yT); yt1(:, 2:2:end, :) = cmpl2(yT);
%     yf1t2 = fftshift(fft(yt1, nt1, 1), 1);
%     yf1t2 = real(yf1t2(:, 1:2:end, :)) + 1i*real(yf1t2(:, 2:2:end, :));
%     yT = ifft(ifftshift(yf1t2, 1), nt1, 1);     % yTc = fftshift(yTc, 1);        % Move the sample corresponding to t=0 to the center of array
%     t1 = ([-floor(nt1/2):floor(nt1/2)])'*dt1; t1 = t1(1:end-1);       % Time in seconds
%     t1 = ifftshift(t1);
%     
%     [tx, ty] = meshgrid(t2, t1);
%     [fx, fy] = meshgrid(del2, del1);
%     %% Construct the model signals for plotting
%     n_peaks = size(true_chsh, 1);    % Total number of peaks
%     omega = 2*pi*bsxfun(@minus, bsxfun(@times, true_chsh', [c1_ref; c2_ref]), [f1_ref; f2_ref]);
%     tau = 0; theta = 0;
%     Z = NaN(nt1, nt2, n_peaks);
%     ZZ = NaN(n_peaks, n_peaks);
%     for i = 1 : n_peaks
%         Z(:,:,i) = kron(exp(1i*omega(1, i)*t1 - true_alpha(i, 1)*abs(t1)), exp(1i*omega(2, i)*t2' - true_alpha(i, 2)*abs(t2')));
%         Z(:,:,i) = Z(:,:,i)/norm(Z(:,:,i));
%         for j = 1 : i
%             ZZ(i, j) = reshape(Z(:,:,i), [], 1)' * reshape(Z(:,:,j), [], 1);
%             ZZ(j, i) = ZZ(i, j)';
%         end;
%     end;
%     ZZ = (ZZ + ZZ')/2;   % Make sure ZZ is Hermitian to avoid numerical errors
%     Zy = reshape(Z, [], n_peaks)' * reshape(yT(:, :, 2), [], 1);
%     b = inv(ZZ)*Zy; b = b / norm(b);     % Inferred b
%     Zb = bsxfun(@times, Z, reshape(b, 1, 1, []));
%     xT = sum(Zb, 3);
%     %% Compute spectra for plotting
%     nt1half = ceil(numel(t1)/2);    % Number of the first-dimension samples in the first half
%     yT_zf = [yT(1:nt1half, :, :); zeros(nf1-nt1, nt2, n_echo); yT(nt1half+1:end, :, :)];   % Zero-filling
%     xT_zf = [xT(1:nt1half, :); zeros(nf1-nt1, nt2); xT(nt1half+1:end, :)];   % Zero-filling
%     % Compute the FT
%     yF = fftshift(fftshift(fft2(yT_zf, nf1, nf2), 1), 2) / sqrt(nt1*nt2);      % frequency-domain representation of y
%     xF = fftshift(fftshift(fft2(xT_zf, nf1, nf2), 1), 2) / sqrt(nt1*nt2);      % frequency-domain representation of y
%     peaks_plot_1 = real(yF(:,:,indY));
%     peaks_plot_2 = real(xF);
% elseif strcmp(mode, 'bicomplex')
%     %% Construct the model signal
%     n_peaks = size(true_chsh, 1);    % Total number of peaks
%     omega = 2*pi*bsxfun(@minus, bsxfun(@times, true_chsh', [c1_ref; c2_ref]), [f1_ref; f2_ref]);
%     Z = bicomplex(NaN(nt1*nt2, n_peaks), NaN(nt1*nt2, n_peaks));
%     ZZ = bicomplex(NaN(n_peaks, n_peaks), NaN(n_peaks, n_peaks));
%     for i = 1 : n_peaks
%         S1 = exp(1i*omega(1, i)*t1 - true_alpha(i, 1)*abs(t1));
%         S2 = exp(1i*omega(2, i)*t2' - true_alpha(i, 2)*abs(t2'));
%         S1z = kron(real(S1), real(S2))+1i*kron(imag(S1), real(S2));
%         S2z = kron(real(S1), imag(S2))+1i*kron(imag(S1), imag(S2));
%         Z(:,i) = bicomplex(reshape(S1z, [], 1), reshape(S2z, [], 1));
%         Z(:,i) = Z(:,i) / norm(norm(Z(:,i)));
%         for j = 1 : i
%             ZZ(i, j) = sum(conj(Z(:,i)) .* Z(:,j));
%             ZZ(j, i) = conj(ZZ(i, j));
%         end;
%     end;
%     ZZ = (ZZ + transpose(conj(ZZ)))/2;   % Make sure ZZ is Hermitian to avoid numerical errors
%     Zy = transpose(conj(Z)) * reshape(yT(:, :, 1), [], 1);
%     ZZi = inv(ZZ);
%     b = ZZi*Zy;     % Inferred b
%     Zb = Z * b;
%     xT = reshape(Zb, nt1, nt2);                     % Modeled signal
%     %% Compute spectra for plotting
% %     % Absolute values of bicomplex data
% %     peaks_plot_1 = abs(modc(ffthc(yT(:,:,indY), [nf1, nf2], 'hann', 'bicomplex')));
% %     peaks_plot_2 = abs(modc(ffthc(xT(:,:), [nf1, nf2], 'hann', 'bicomplex')));
%     
%     % Real spectra of reduced complex data
%     peaks_plot_1 = real(ffthc(yT(:,:,indY), [nf1, nf2], 'rcos', 'complex'));
%     peaks_plot_2 = real(ffthc(xT(:,:), [nf1, nf2], 'rcos', 'complex'));
% end;
%% Compute the SNR
% nT = bicomplex(sig_est*(randn(nt1, nt2, n_echo) + 1i*randn(nt1, nt2, n_echo)),...
%     sig_est*(randn(nt1, nt2, n_echo) + 1i*randn(nt1, nt2, n_echo)));
% snrEst = 10*log10(real(sum(sum(sum(conj(xT(:,:,1)).*xT(:,:,1), 1), 2), 3)) / real(sum(sum(sum(conj(nT(:,:,1)).*nT(:,:,1), 1), 2), 3)) );     % snrAllTrue(ii, p, q, j) = 10*log10(real(transpose(conj(xT(:)))*xT(:)) / real(transpose(conj(nT(:)))*nT(:)));

%% Plot the spectra
indx_plot = 1:n_peaks;    % Which peaks to mark and label
indx_label = find(and(orig_chsh(:,1)>7.8, orig_chsh(:,1)<8.5) .* and(orig_chsh(:,2)>104, orig_chsh(:,2)<129));

figure('Position', [180, 120, 1050, 530]);
% set(gcf, 'Position', [450 115 745 850]);
ax(1) = subplot(1, 2, 1);
contour(del1, del2, fliplr(rot90(peaks_plot_1, -1)), linspace(0.75*min(min(peaks_plot_1)), 0.75*max(max(peaks_plot_1)), 100));    % 0.5*min(min(peaks_plot_1))
hold all;
scatter(true_chsh(indx_plot,1), true_chsh(indx_plot,2), 3, 'r*');
text(true_chsh(indx_label, 1)-0.0001, true_chsh(indx_label, 2)-0.001*3.5, num2str(indx_label(:)), 'FontSize', 7);
set(gca,'XDir','reverse'); set(gca,'YDir','reverse');
ylabel('^1^5N, ppm'); xlabel('^1H, ppm');
xlim([6.5 10]); ylim([104 129]);
title('Observed spectrum');

ax(2) = subplot(1, 2, 2);
contour(del1, del2, fliplr(rot90(peaks_plot_2, -1)), linspace(0.5*min(min(peaks_plot_2)), 0.65*max(max(peaks_plot_2)), 42));
hold all;
line([true_chsh(indx_plot,1) best_chsh(indx_plot,1)]', [true_chsh(indx_plot,2) best_chsh(indx_plot,2)]');
scatter(true_chsh(indx_plot,1), true_chsh(indx_plot,2), 7, 'ro');
scatter(best_chsh(indx_plot,1), best_chsh(indx_plot,2), 3, 'r*');
text(best_chsh(indx_label, 1)-0.0001, best_chsh(indx_label, 2)-0.001*3.5, num2str(indx_label(:)), 'FontSize', 7);
% for i = 1 : n_peaks
% rectangle('Position', [lwr_pars(i, 1), lwr_pars(i,2), upr_pars(i,1)-lwr_pars(i,1), upr_pars(i,2)-lwr_pars(i,2)]);
% end;
set(gca,'XDir','reverse'); set(gca,'YDir','reverse');
ylabel('^1^5N, ppm'); xlabel('^1H, ppm');
xlim([6.5 10]); ylim([104 129]);
title('Modeled spectrum');

linkaxes(ax);
return;

%% Plot a 2D objective function
clear logPostHCph;
k1 = 1+n_peaks*(4-1); xlbl = 'Relaxation rate';    % Which dimension to plot
k2 = 15+n_peaks*(4-1); ylbl = 'Relaxation rate';    % Which dimension to plot
x = linspace(lwr_pars(k1), upr_pars(k1), 100);
y = linspace(lwr_pars(k2), upr_pars(k2), 100);
[xx, yy] = meshgrid(x, y);
for i = 1 : numel(x)
    for j = 1 : numel(y)
        fprintf('Computing for i = %d, j = %d.\n', i, j);
        pars = allPtcl(:, end);
        pars(k1) = x(i);
        pars(k2) = y(j);
        P_smpl(i,j) = logPostHCph(yT1, pars, c_ref, f_ref, t1, t2, 'tau', allTau(:, end), 'noPhasing', 0);
    end;
end;

figure('Position', [300 320 880 330]);
contour(x, x, P_smpl, 25);
xlabel(xlbl); ylabel(ylbl);

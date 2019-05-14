%%
% This script provides an example of parameters fitting. Please note that
% fitting spectra with multiple peaks may take significant time, and manual
% inspection of the fitting results is advised.

%% Default parameters
clear;
s = RandStream('mcg16807','Seed',0);
RandStream.setGlobalStream(s);

%% Choose an experiment to run
example = 1150;
switch example
    case 1150  % T1
        load(fullfile('_data', 'fitted_pars_1150_108_ph'));
        t_echo = [0.08 0.16 0.24 0.4 0.56 0.64 0.72 0.8];
    case 1151  % T1rho
        load(fullfile('_data', 'fitted_pars_1151_108_ph'));
        t_echo = [0.001 0.005 0.01 0.018 0.035 0.055 0.08 0.11 0.15 0.2];
    case 1155  % T2
        load(fullfile('_data', 'fitted_pars_1155_108_ph'));
        t_echo = [0.0144 0.0288 0.0432 0.0576 0.072 0.0864 0.1008 0.1152 0.1296 0.144];
end;

% Load a specific Azara file (1150-T1, 1151-T1rho, 1155-T2)
[yT, c_ref, f_ref, t, del] = loadData_115x(example);

% Initialize the results
best_chsh = orig_chsh; best_alpha = orig_alpha;

% Only analyze the first 104 peaks
best_chsh = best_chsh(1:104, :);
best_alpha = best_alpha(1:104, :);

nt = size(yT); nt1 = nt(1); nt2 = nt(2);      % Number of points in each time dimension
n_peaks = size(best_chsh, 1);    % Total number of peaks
n_echo = numel(t_echo);          % Number of relaxation planes

%% Fit the model to the data
%% Step 1. Initialization of the parameters
best_pars = [best_chsh best_alpha];

%% Step 2. Set the range limits for optimization
D = pdist2(best_chsh, best_chsh); Dmin = D; Dmin(1:n_peaks+1:end) = Inf;
Dmin = min(Dmin)'; Dmin = min(Dmin, 0.5);   % Distance to the closest peak
lwr_pars = bsxfun(@minus, best_pars, [min(Dmin, 0.05)/2 min(Dmin, 0.5)/2 35*ones(n_peaks,1) 30*ones(n_peaks,1)]);
lwr_pars(:, 3) = 0; lwr_pars(:, 4) = 0;
upr_pars = bsxfun(@plus, best_pars, [min(Dmin, 0.05)/2 min(Dmin, 0.5)/2 20*ones(n_peaks,1) 20*ones(n_peaks,1)]);
upr_pars(:, 3) = 200; upr_pars(:, 4) = 200;

%% Step 3. Define the indexes of parameters to be optimized. The indices
% will be grouped in the cell array with each cell containing indices of
% parameters optimized simultaneously.
% Indices in the ranges 1...n_peaks and n_peaks+1...2*n_peaks corrspond to
% chemical shifts in both dimensions.
% Indices 2*n_peaks+1...3*n_peaks and 3*n_peaks+1...4*n_peaks correspond to
% relaxation rates alpha.
indx_optm = reshape(1:(n_peaks*4), [], 4);

% For example, fit only chemical shifts of each peak separately
indx_optm = indx_optm(:, 1:2);
indx_optm = indx_optm';
indx_optm = indx_optm(:)';
indx_optm = mat2cell(indx_optm, 1, 2*ones(1, n_peaks));
indx_optm = indx_optm([randperm(numel(indx_optm)) randperm(numel(indx_optm)) randperm(numel(indx_optm))]);

% % % Or, as another example, fit all parameters of each peak 
% % indx_optm = indx_optm([203], 1:4);
% % indx_optm = indx_optm';
% % indx_optm = indx_optm(:)';
% % indx_optm = mat2cell(indx_optm, 1, 4*ones(1, numel(indx_optm)/4));

%% Step 4. Loop over (some or all) cells in the indx_optm array
numTryPts = 0;    % Number of trial points for MATLAB global search. Set >0 to try several starting points for optimization

for kk = [4, 9, 10]      %    [18    20    39    45    48    49    71    73   100   119   121   123   152   153   154   172   175   207] % randperm(n_peaks)    %  1:numel(indx_optm)
    fprintf('Analyzing step #%d.\n', kk);
    pars_new = best_pars(:)';

% ----- Form parameter-selecting matrices
    K = ones(1, 4*n_peaks); K(indx_optm{kk}) = 0;
    K = diag(K);   % Keep
    if nnz(K) == 0; K = 0; O = 1;
    else
        O = eye(4*n_peaks) - K; O = O(indx_optm{kk}, :);   % Optimize
        K = pars_new*K;
    end;
    costFunc_rdcd = @(x) costFunc_VP(yT(:,:,1), x*O+K, c_ref, f_ref, t{1}, t{2});
    pars_0 = [pars_new(indx_optm{kk})];
   
% --- MATLAB global search
    opts = optimoptions(@fmincon,'Algorithm','sqp');
    if numTryPts == 0
        [gsPtcl, gsPost, flg, output] = fmincon(@(x) costFunc_rdcd(x), pars_0, ...
            [],[],[],[],lwr_pars(indx_optm{kk}),upr_pars(indx_optm{kk}), [], opts);
    else
        gsProblem = createOptimProblem('fmincon','objective',@(x) costFunc_rdcd(x),...
            'x0', pars_0, 'lb', [lwr_pars(indx_optm{kk})], 'ub', [upr_pars(indx_optm{kk})], 'options', opts);
        gs = GlobalSearch('Display', 'iter', 'NumTrialPoints', 10, 'NumStageOnePoints', 5, 'DistanceThresholdFactor', 0.9);    % gs = GlobalSearch('Display', 'iter', 'NumTrialPoints', 250, 'NumStageOnePoints', 150, 'DistanceThresholdFactor', 0.9);
        [gsPtcl, gsPost, flg, og, solns] = run(gs,gsProblem);
    end;
    pars_fit = gsPtcl*O + K;
    
    best_alpha = pars_fit(2*n_peaks+1:4*n_peaks);
    best_alpha = reshape(best_alpha, [], 2);
    best_chsh = pars_fit(1:2*n_peaks);
    best_chsh = reshape(best_chsh, [], 2);
    fprintf('\n');
end;

%% Step 5. Found best parameters
best_pars = reshape(pars_fit, [], 4);


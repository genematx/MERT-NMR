%%
% This script presents examples of analysing 2D bicomplex NMR data in the
% time domain. Please select a specific example below.
% Required inputs are:
% 1. Time-domain relaxationdata in the Azara format
% 2. Matrices of model parameters (chemical shifts and peak widths)
%    estimated, for example during peak picking. Additionally, acquisition
%    delays in both dimensions can be supplied for propper phasing of the
%    spectrum for plotting.
% 3. An array of relaxation times t_echo.
% 
% The result of computations is the array of relation rates along with
% plots of modelled and experimntal spectra.

%% Default parameters
clear;
s = RandStream('mcg16807','Seed',0);
RandStream.setGlobalStream(s);
addpath('_tools\whittaker_smoother');   % Whittaker smoother for baseline calculation

nf1 = 1024;        nf2 = 1024;       % Number of points in each dimension of the spectrum (after zero-filling)
indY = 1;             % Use as the starting plane to plot
xlims = [6.5 10]; ylims = [104 129];    % Limits for plotting the spectra

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

% % Uncomment to use the original parameters determined by peak-picking without optimization
% best_chsh(1:size(orig_chsh,1), :) = orig_chsh;
% best_alpha(1:size(orig_alpha, 1), :) = orig_alpha;

% % Only analyze the first 104 peaks
% best_chsh = best_chsh(1:104, :);
% best_alpha = best_alpha(1:104, :);

nt = size(yT); nt1 = nt(1); nt2 = nt(2);      % Number of points in each time dimension
n_peaks = size(best_chsh, 1);    % Total number of peaks
n_echo = numel(t_echo);          % Number of relaxation planes

%% Fit the baseline
% Construct the model signal
Z = getModelMatrix(best_chsh, best_alpha, c_ref, f_ref, t{1}, t{2}, 'tau', best_tau);
ZZ = transpose(conj(Z))*Z;
Zy = transpose(conj(Z)) * reshape(yT(:, :, indY), [], 1);
ZZi = inv(ZZ);
b = ZZi*Zy;     % Inferred vector of intensities b

xT = reshape(Z * b, nt1, nt2);                     % Ideal modeled signal in the time domain
yF = ffthc(yT(:,:,indY), [nt1, nt2], 'none', 'bicomplex');     % Measured spectrum
xF = ffthc(xT, [nt1, nt2], 'none', 'bicomplex');               % Modelled spectrum
rF = yF - xF;                                                  % Residual spectrum

% Compute the baseline by smoothing the residual
% % Case 1. Smoothen in both dimensions
% bF = bicomplex((wsmooth(real(rF)) + 1i*wsmooth(imag1(rF))), ...
%     (wsmooth(imag2(rF)) + 1i*wsmooth(imag12(rF))) );
% Case 2. Smoothen only in the direct dimension
bF = bicomplex(zeros(nt1, nt2), zeros(nt1, nt2) );
for i = 1 : nt2
    bF(:, i) = bicomplex((wsmooth(real(rF(:, i))) + 1i*wsmooth(imag1(rF(:, i)))), ...
    (wsmooth(imag2(rF(:, i))) + 1i*wsmooth(imag12(rF(:, i)))) );
end;
bT = iffthc(bF, [nt1, nt2]);
% % Case 3. Constant baseline
% bT = bicomplex([1; zeros(nt1*nt2-1, 1)], zeros(nt1*nt2, 1));

% Add the baseline to the model matrix, renormalize and update ZZ
Z = [Z bT(:)];
Z = Z*diag(1./sqrt(diag(real(transpose(conj(Z))*Z))));     % Normalize Z
ZZ = transpose(conj(Z))*Z;
ZZi = inv(ZZ);

%% Construct the model signal
ZY = transpose(conj(Z)) * reshape(yT, [], n_echo);
b = ZZi*ZY;     % Inferred b
Zb = Z * b;                  % Plot WITH baseline
%    Zb = Z(:, 1:n_peaks) * b(1:n_peaks,:);             % Plot WITHOUT baseline
xT = reshape(Zb, nt1, nt2, n_echo);                     % Modeled signal

%% Compute spectra for plotting
% % Absolute values of bicomplex data
%peaks_plot_1 = abs(modc(ffthc(yT(:,:,:), [nf1, nf2], 'hann', 'bicomplex')));
%peaks_plot_2 = abs(modc(ffthc(xT(:,:,:), [nf1, nf2], 'hann', 'bicomplex')));

% Real spectra of reduced complex data
peaks_plot_1 = real(ffthc(yT(:,:,:), [nf1, nf2], 'sinebell', 'complex'));
peaks_plot_2 = real(ffthc(xT(:,:,:), [nf1, nf2], 'sinebell', 'complex'));
%% Plot the spectra
n_peaks = size(best_chsh, 1);    % Total number of peaks
indx_plot = 1:n_peaks;    % Which peaks to mark and label
indx_label = 1:n_peaks;
del1 = del{1}; del2 = del{2};

lvl_min = 0.04*max([peaks_plot_1(:); peaks_plot_2(:)]);
lvl_max = 0.85*max([peaks_plot_1(:); peaks_plot_2(:)]);

figure('Position', [180, 120, 1050, 530]);
ax(1) = subplot(1, 2, 1);
grid on; hold all;
contour(del1, del2, fliplr(rot90(peaks_plot_1(:,:,indY), -1)), linspace(lvl_min, lvl_max, 8));
scatter(best_chsh(indx_plot,1), best_chsh(indx_plot,2), 3, 'r*');
text(best_chsh(indx_label, 1)+0.001*10, best_chsh(indx_label, 2)+0.001*3.5, num2str(indx_label(:)), 'FontSize', 7);
set(gca,'XDir','reverse'); set(gca,'YDir','reverse');
ylabel('^1^5N, ppm'); xlabel('^1H, ppm');
xlim(xlims); ylim(ylims);

title('Observed spectrum');

ax(2) = subplot(1, 2, 2);
grid on; hold all;
contour(del1, del2, fliplr(rot90(peaks_plot_2(:,:,indY), -1)), linspace(lvl_min, lvl_max, 8));
line([best_chsh(indx_plot,1) best_chsh(indx_plot,1)]', [best_chsh(indx_plot,2) best_chsh(indx_plot,2)]');
scatter(best_chsh(indx_plot,1), best_chsh(indx_plot,2), 5, 'ro');
scatter(best_chsh(indx_plot,1), best_chsh(indx_plot,2), 3, 'r*');
text(best_chsh(indx_label, 1)+0.001*10, best_chsh(indx_label, 2)+0.001*3.5, num2str(indx_label(:)), 'FontSize', 7);
set(gca,'XDir','reverse'); set(gca,'YDir','reverse');
ylabel('^1^5N, ppm'); xlabel('^1H, ppm');
xlim(xlims); ylim(ylims);
title('Modeled spectrum');

linkaxes(ax);
% return;

%% Plot a single peak in 1D
indPlt_peak = 37;            % Choose which peak to plot
indPlt_plane = [1 3 5 7];    % Choose which planes to show
% Find the indices in the plotting arrays that correspond to the chosen peak
[~, ind1] = min(abs(del1 - best_chsh(indPlt_peak,1)));
[~, ind2] = min(abs(del2 - best_chsh(indPlt_peak,2)));

figure('Position', [500, 100, 900, 850]);
title(['Peak #', int2str(indPlt_peak)]);
for i = 1 : 4
    pl1(i) = subplot(4,2,(i-1)*2+1);
    plot(del2, peaks_plot_1(ind1,:,indPlt_plane(i)));
    hold all
    plot(del2, peaks_plot_2(ind1,:,indPlt_plane(i)));
    xlim(ylims);
    ylim([-5, 65]);
    plot(best_chsh(indPlt_peak,2), 0, 'g*', 'MarkerSize', 1);
    ylabel(['Plane #', int2str(indPlt_plane(i))]);
    if i == 1
        legend('Measured data', 'Fitted model')
    end;
    
    pl2(i) = subplot(4,2,i*2);
    plot(del1, peaks_plot_1(:,ind2,indPlt_plane(i)));
    hold all
    plot(del1, peaks_plot_2(:,ind2,indPlt_plane(i)));
    xlim(xlims);
    ylim([-5, 65]);
    plot(best_chsh(indPlt_peak,1), 0, 'g*', 'MarkerSize', 1);
end;
linkaxes(pl1); linkaxes(pl2);

% % Plot residuals
% figure('Position', [500, 100, 900, 850]);
% title(['Peak #', int2str(indPlt_peak)]);
% for i = 1 : 4
%     pl1(i) = subplot(4,2,(i-1)*2+1);
%     plot(del2, peaks_plot_1(ind1,:,indPlt_plane(i)) - peaks_plot_2(ind1,:,indPlt_plane(i)));
%     hold all
%     xlim(ylims);
%     ylim([-4.5, 4.5]);
%     plot(best_chsh(indPlt_peak,2), 0, 'g*', 'MarkerSize', 1);
%     ylabel(['Plane #', int2str(indPlt_plane(i))]);
%     
%     pl2(i) = subplot(4,2,i*2);
%     plot(del1, peaks_plot_1(:,ind2,indPlt_plane(i)) - peaks_plot_2(:,ind2,indPlt_plane(i)));
%     hold all
%     xlim(xlims);
%     ylim([-4.5, 4.5]);
%     plot(best_chsh(indPlt_peak,1), 0, 'g*', 'MarkerSize', 1);
% end;
% linkaxes(pl1); linkaxes(pl2);
% %return;

%% Estimate T1/T2/T1rho
% Estimate peak intensities
w = ZZi * ZY;
w = sqrt(real(conj(w).*w));

indx_planes = [1:n_echo];    % Which planes to consider for fitting
for i = 1 : n_peaks
    % Method 1. Fit a linear function to the logarithms
    linFit = polyfit(t_echo(indx_planes), log(w(i, indx_planes)), 1);
    beta_hat_linFit(i) = -linFit(1);
    % Method 2. Fit an exponential directly
    expFit = fit(t_echo(indx_planes)',w(i, indx_planes)','exp1', 'Exclude', isnan(w(i, indx_planes)));
    beta_hat_expFit(i) = -expFit.b;
end;
beta_hat = beta_hat_expFit(:);

% Comparison figure
R = corrcoef(orig_beta(~isnan(orig_beta)), beta_hat(~isnan(orig_beta)));
figure;
%scatter(orig_beta, beta_hat, 'b.');
plot(orig_beta(1:numel(beta_hat)), beta_hat, 'r^', 'MarkerSize', 3);
hold all;
% scatter(orig_beta, beta_hat_expFit(:), 'b.');
line([0 15], [0 15]);
hold all;
xlabel('FT integration');
ylabel('Model-based');
axis square;
grid on;
% xlim([0.5 3.0]); ylim([0.5 3.0]);    % For T1
% xlim([0 10]); ylim([0 10]);          % For T1rho
% xlim([1 11]); ylim([1 11]);          % For T2

function result = ffthc(yT, NF, wnd, sptype)
% Computes the FT of a bicomplex-valued matrix
% yT - time-domain data, 2D matric of cmplex or bicomplex values
% NF = [nf1, nf2] - sizes of the resulting FT
% wnd - a string specifying a windowing function to be applied; can be
%       'none', 'gauss', 'tukey', 'hann', 'rcos', 'sinebell', 'exp', 'qsine'
% sptype - type of the resulting spectrum (either 'complex' or 'bicomplex')
% 
% phase = [ph0 ph1] - zero and first order phase correction terms applied
%       along the direct dimension

NT = size(yT);
if size(NT, 2) > 2
    result = NaN(NF(1), NF(2), NT(3));
    if strcmp(sptype, 'bicomplex'); result = bicomplex(result, result); end;
    for i = 1 : NT(3)
        result(:,:,i) = ffthc(yT(:,:, i), NF, wnd, sptype);
    end;
    return;
%     error('The input must be two-dimensional.');
end;
nt1 = NT(1); nt2 = NT(2);
nf1 = NF(1); nf2 = NF(2);

% Define weighting windows in the time domains
switch wnd
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
    case 'rcos'      % Raised cosine window
        t = linspace(0, pi+pi/3, nt1)'; w1 = cos(t - pi/3)+1;
        t = linspace(0, pi+pi/3, nt2)'; w2 = cos(t - pi/3)+1;
    case 'sinebell'  % Sinebell function with 60 deg offset
        w1 = cos(linspace(-pi/3, pi/2, nt1)');
        w2 = cos(linspace(-pi/3, pi/2, nt2)');
    case 'exp'      % Exponential window
        t = linspace(0, 3*exp(1), nt1)'; w1 = exp(-1*t);
        t = linspace(0, 3*exp(1), nt2)'; w2 = exp(-1*t);
    case 'qsine'     % Raised and squared cosine window
        t = linspace(0, pi+pi/3, nt1)'; w1 = ((cos(t - pi/3)+1)/2).^2;
        t = linspace(0, pi+pi/3, nt2)'; w2 = ((cos(t - pi/3)+1)/2).^2;
    otherwise
        error('Unknown window name.');
end;

% Rearrange data as a complex matrix with alternating R/I columns
yt1 = NaN(nt1, 2*nt2);     % Data in time domain, complex values along the 1st dimenasion, R-I-R-I-... along the second
yt1(:, 1:2:end) = cmpl1(yT(:,:));
yt1(:, 2:2:end) = cmpl2(yT(:,:));

% Compute the FT in the first dimension
yf1t2 = fftshift(fft(bsxfun(@times, yt1, w1), nf1, 1), 1);
yt2r = real(yf1t2(:, 1:2:end)) + 1i*real(yf1t2(:, 2:2:end));
yt2i = imag(yf1t2(:, 1:2:end)) + 1i*imag(yf1t2(:, 2:2:end));
% Compute the FT in the sewcond dimension
yf2r = fftshift(fft(bsxfun(@times, yt2r, w2'), nf2, 2), 2);
yf2i = fftshift(fft(bsxfun(@times, yt2i, w2'), nf2, 2), 2);

switch sptype
    case 'complex'    % Phased complex spectrum
        result = yf2r;
    case 'bicomplex'
        result = bicomplex(yf2r, yf2i);
    otherwise
        error('Unknown type of the spectrum.');
end;

end


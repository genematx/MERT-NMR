function result = iffthc(yF, NT)
% Computes the iFFT of a bicomplex-valued matrix
% yF - frequency-domain data, 2D matric of complex or bicomplex values
% NT = [nt1, nt2] - sizes of the resulting iFT

NF = size(yF);
if size(NF, 2) > 2; error('The input must be two-dimensional.'); end;
nt1 = NT(1); nt2 = NT(2);
nf1 = NF(1); nf2 = NF(2);

yf2r = cmpl1(yF);
yf2i = cmpl2(yF);
yt2r = ifft(ifftshift(yf2r, 2), nt2, 2);
yt2i = ifft(ifftshift(yf2i, 2), nt2, 2);
yf1t2 = NaN(nf1, 2*nt2);
yf1t2(:, 1:2:end) = real(yt2r) + 1i*real(yt2i);
yf1t2(:, 2:2:end) = imag(yt2r) + 1i*imag(yt2i);
yt1 = ifft(ifftshift(yf1t2, 1), nt1, 1);
result = bicomplex(yt1(:, 1:2:end), yt1(:, 2:2:end));

end


function zeta = mat2bicomp(mat)
% Takes the matrix representation and returns the corresponding bicomplex
% structure
sizes = size(mat);
str1 = '1:sizes(1)/2,1:sizes(2)/2';
str2 = 'sizes(1)/2+1:end,1:sizes(2)/2';
for i = 3:length(sizes);
    str1 = [str1 sprintf(',1:sizes(%i)',i)];
    str2 = [str2 sprintf(',1:sizes(%i)',i)];
end
str1 = sprintf('mat(%s)',str1);
str2 = sprintf('mat(%s)',str2);
z1 = eval(str1);
z2 = eval(str2);
zeta = bicomplex(z1, z2);
end
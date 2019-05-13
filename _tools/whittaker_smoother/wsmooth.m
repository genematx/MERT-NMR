function [zs,L] = wsmooth(z,x,y,L)
%WSMOOTH   1D and 2D robust smoothing.
%   ZS = WSMOOTH(Z,X,L) smoothes the signal Z(X) using the smoothing
%   parameter L. Z and X must be vector arrays of same length. L must be a
%   negative or positive real scalar. The larger L is, the smoother the
%   output will be. Typical values for L are between 0 and 3.
%
%   ZS = WSMOOTH(Z,X,Y,L) smoothes the surface Z(X,Y). X and Y must be
%   vector arrays whose lengths match the size of matrix Z.
%
%   ZS = WSMOOTH(Z,L) uses unit spacings for X and/or Y.
%
%   ZS = WSMOOTH(Z,X,Y) or ZS = WSMOOTH(Z,X) or ZS = WSMOOTH(Z): if L is
%   omitted, it is automatically determined using an incremental process
%   based on an estimation of the output signal variance. The function
%   ESTIMATENOISE written by John D'Errico is required for this option
%   (http://www.mathworks.com/matlabcentral/files/16683/estimatenoise.m).
%   [ZS,L] = WSMOOTH(...) also returns the proposed value for L so that you
%   can fine-tune the smoothing subsequently if needed.
%   
%   Notes:
%   -----
%   Smoothing also works if some data are missing (see examples); data are
%   considered missing if they are not finite (NaN or Inf). WSMOOTH can
%   thus be used as an alternative to INPAINT_NANS (another John D'Errico's
%   function): Z = INPAINT_NANS(ZNAN) and Z = WSMOOTH(ZNAN,L) with L<=0
%   give nearly similar results. Smoothing level, however, can be adjusted
%   with WSMOOTH if necessary.
%
%   This program is based on the Whittaker's smoother. The 1D algorithm has
%   been well described by Eilers in Anal. Chem. 2003, 75, 3631-3636. This
%   algorithm has been here adapted to 2D. Note that a 3D update could be
%   easily written.
%   
%   Examples:
%   --------
%   % 1-D example
%   x = 1:100;
%   y = cos(x/10)+(x/50).^2;
%   yn = y + 0.2*randn(1,100);
%   ys = wsmooth(yn,2);
%   plot(x,yn,'.-',x,ys)
%
%   % 1-D example with missing data and non uniform grid
%   x = sort(rand(1,100)*99)+1;
%   y = cos(x/10)+(x/50).^2;
%   y(45:60) = NaN;
%   yn = y + 0.2*randn(1,100);
%   ys = wsmooth(yn,x);
%   plot(x,yn,'.-',x,ys)
%
%   % 2-D example
%   [xp,yp] = deal(0:.02:1);
%   [x,y] = meshgrid(xp,yp);
%   f = exp(x+y) + sin((x-2*y)*3);
%   fn = f + (rand(size(f))-0.5);
%   fs = wsmooth(fn);
%   subplot(121), surf(xp,yp,fn), zlim([0 8])
%   subplot(122), surf(xp,yp,fs), zlim([0 8])
%
%   % 2-D example with missing data and non-uniform grid
%   xp = [0 sort(rand(1,48)) 1]; yp = [0 sort(rand(1,48)) 1];
%   [x,y] = meshgrid(xp,yp);
%   f = exp(x+y) + sin((x-2*y)*3);
%   f(round(rand(1,100)*2500)) = NaN; f(20:30,20:40) = NaN;
%   fn = f + (rand(size(f))-0.5);
%   fs = wsmooth(fn,xp,yp,1.5);
%   subplot(121), surf(xp,yp,fn), zlim([0 8])
%   subplot(122), surf(xp,yp,fs), zlim([0 8])
%
%   % comparison with inpaint_nans when removing NaN
%   [x,y] = meshgrid(0:.02:1);
%   [znan,z0] = deal(exp(x+y));
%   znan(10:25,20:35) = NaN;
%   znan(15:45,2:5) = NaN;
%   znan(35:37,20:45) = NaN;
%   z1 = inpaint_nans(znan);
%   z2 = wsmooth(znan,0);
%   subplot(121), surf(znan), zlim([1 8])
%   subplot(122), surf(z2), zlim([1 8]), title('wsmooth')
%
%   See also SMOOTH, ESTIMATENOISE, INPAINT_NANS.
%
%   -- Damien Garcia -- 2008/03
%   -- Special thanks to John D'Errico (for ESTIMATENOISE)


% % Check input arguments
% ---
narginchk(1,4)

 if ndims(z)~=2
     error('MATLAB:wsmooth:WrongZArray',...
         'First input must be a vector or a matrix.')
 end

if nargin==1 % wsmooth(z)
    if isvector(z) % 1D
        z = z(:)';
        nx = length(z);
        x = 1:nx; y = 1;
    elseif ndims(z)==2 % 2D
        [ny,nx] = size(z);
        x = 1:nx; y = 1:ny;
    else
        error('MATLAB:wsmooth:WrongInputs',...
            'Wrong input arguments: type ''help wsmooth''.')
    end
elseif nargin==2
    if isscalar(x) && ~isvector(z) % wsmooth(z,L), 2D
        L = x;
        [ny,nx] = size(z);
        x = 1:nx; y = 1:ny;
    elseif isscalar(x) && isvector(x) % wsmooth(z,L), 1D
        L = x;
        z = z(:)';
        nx = length(z);
        x = 1:nx; y = 1;
    elseif isvector(x) && isvector(z) % wsmooth(z,x), 1D
        z = z(:)';
        y = 1;
    else
        error('MATLAB:wsmooth:WrongInputs',...
            'Wrong input arguments: type ''help wsmooth''.')
    end
elseif nargin==3
    if ~isvector(z) && isvector(x) && isvector(y) % wsmooth(z,x,y), 2D
        % Nothing to declare
    elseif isvector(z) && isvector(x) && isscalar(y) % wsmooth(z,x,L), 1D
        L = y;
        z = z(:)';
        y = 1;
    else
        error('MATLAB:wsmooth:WrongInputs',...
            'Wrong input arguments: type ''help wsmooth''.')
    end
else % wsmooth(z,x,y,L), 2D
    if isvector(z) || ~isvector(x) || ~isvector(y) || ~isscalar(L)
        error('MATLAB:wsmooth:WrongInputs',...
            'Wrong input arguments: type ''help wsmooth''.')   
    end
end

[ny,nx] = size(z);
if ~isequal(length(x),nx) || ~isequal(length(y),ny)
    error('MATLAB:wsmooth:XYLengthMismatch',...
        'Number of elements in X and/or Y is inadequate.')
elseif any(sign(diff(x))~=1) && any(sign(diff(x))~=-1)
    error('MATLAB:wsmooth:NonMonotonicX',...
        'X must be strictly monotonic.')
elseif ~isequal(y,1) && any(sign(diff(y))~=1) && any(sign(diff(y))~=-1)
    error('MATLAB:wsmooth:NonMonotonicY',...
        'Y must be strictly monotonic.')
end

if exist('L','var')
    L = 10^L;
elseif ~exist('estimatenoise','file')
    error('MATLAB:wsmooth:MissingFunction',...
        ['ESTIMATENOISE is required if the L parameter is not specified.'...
        '\nESTIMATENOISE can be downloaded from the Matlab FEX (file ID: 16683).'])
end

x = x(:);
y = y(:);
class0 = class(z);
z = double(z);

% % 1D
if isequal(y,1)
    % Create the X-difference matrix
    % FDM: 2nd derivative, X-central difference, 2nd order
    dx0 = x(2:end-1)-x(1:end-2);
    dx1 = x(3:end)-x(2:end-1);    
    D = spdiags([2./dx0./(dx0+dx1) -2./dx0./dx1 2./dx1./(dx0+dx1)],...
        [0 1 2],nx-2,nx)*mean(dx0.^2);
    
    % Weight function (0 for NaN/Inf data, 1 otherwise)
    I = isfinite(z);
    z(~I) = 0;
    W = spdiags(double(I(:)),0,nx,nx);
    
    % Solve the linear system
    if exist('L','var')
%         L = sparse(diag(linspace(1, L, nx)));
        zs = (W + L*D'*D)\z(:);
    else
        amp = max(z(I))-min(z(I));
        noisevar = Inf;
        L = -1/3;
        while sqrt(noisevar)/amp > 1e-4
            L = L+1/3;
            zs = (W + 10^L*D'*D)\z(:);
            noisevar = estimatenoise(zs,x);
        end
    end
    zs = reshape(zs,size(z));

% % 2D
else
    % -----
    % The Whittaker's method is adapted to 2D. We solve:
    %   (I + lambda(Dx' Dx + Dy' Dy))zs = z,
    % where Dx and Dy are the difference matrices related to x- and y-
    % direction, respectively.
    % -----
    
    % Create the X-Difference matrix
    % FDM: 2nd derivative, X-central difference, 2nd order
    n = (nx-2)*ny;
    Id = repmat((1:n),[3 1]);
    [Jd,Dx] = deal(zeros(3,n));

    zI = true(ny,nx);
    zI(:,[1 nx]) = false;
    [I,J] = find(zI);

    Jd(1,:) = I+(J-2)*ny;
    Jd(2,:) = I+(J-1)*ny;
    Jd(3,:) = I+J*ny;
        
    dx0 = x(J)-x(J-1);
    dx1 = x(J+1)-x(J);
    Dx(1,:) = 2./dx0./(dx0+dx1);
    Dx(2,:) = -2./dx0./dx1;
    Dx(3,:) = 2./dx1./(dx0+dx1);
    Dx = sparse(Id,Jd,Dx,n,nx*ny)*mean(dx0.^2);

    % Create the Y-Difference matrix
    % FDM: 2nd derivative, Y-central difference, 2nd order
    n = nx*(ny-2);
    Id = repmat((1:n),[3 1]);
    [Jd,Dy] = deal(zeros(3,n));

    zI = true(ny,nx);
    zI([1 ny],:) = false;
    [I,J] = find(zI);

    Jd(1,:) = I-1+(J-1)*ny;
    Jd(2,:) = I+(J-1)*ny;
    Jd(3,:) = I+1+(J-1)*ny;

    dy0 = y(I)-y(I-1);
    dy1 = y(I+1)-y(I);
    Dy(1,:) = 2./dy0./(dy0+dy1);
    Dy(2,:) = -2./dy0./dy1;
    Dy(3,:) = 2./dy1./(dy0+dy1);
    Dy = sparse(Id,Jd,Dy,n,nx*ny)*mean(dy0.^2);
    
    % Weight function (0 for NaN/Inf data, 1 otherwise)
    I = isfinite(z);
    z(~I) = 0;
    W = spdiags(double(I(:)),0,nx*ny,nx*ny);

    % Solve the linear system
    if exist('L','var')
        zs = (W + L*Dx'*Dx + L*Dy'*Dy)\z(:);
        zs = reshape(zs,ny,nx);
    else
        amp = max(z(I))-min(z(I));
        noisevar = Inf;
        L = -1/3;
        while sqrt(noisevar)/amp > 1e-4
            L = L+1/3;
            zs = (W + 10^L*Dx'*Dx + 10^L*Dy'*Dy)\z(:);
            zs = reshape(zs,ny,nx);
            noisevar = estimatenoise(zs,x,2);
        end
    end
end

zs = cast(zs,class0);


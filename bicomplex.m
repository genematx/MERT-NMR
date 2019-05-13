% A class for working with bicomplex numbers in Matlab.
% Created by Adriaen Verheyleweghen, verheyle($at$)stud.ntnu.no
% Please see "Computation of higher-order derivatives using the
% multi-complex step method" for more detail.

classdef bicomplex
    %% BICOMPLEX(z1,z2)
    % Creates an instance of a bicomplex object.
    % zeta = z1 + j*z2, where z1 and z2 are complex numbers.
    
    properties
        z1, z2
    end
    
    methods % Initialization as BC = z1 + j*z2; to initialize the number in polar coordinates, provide three inputs, BC = z1 * z2 * z3 -- z3 can be a scalar.
        function self = bicomplex(varargin)
            switch nargin
                case 2
                    if ~isequal(size(varargin{1}),size(varargin{2}))
                        error('Inputs must be equally sized');
                    end
                    self.z1 = varargin{1};
                    self.z2 = varargin{2};
                case 3
                    if ~isequal(size(varargin{1}),size(varargin{2}))
                        error('The first two inputs must be equally sized');
                    elseif and(~isequal(size(varargin{2}),size(varargin{3})), numel(varargin{3}) ~= 1)
                        error('The size of the third input must match the first two or be a scalar');
                    end
                    self.z1 = real(varargin{1}).*real(varargin{2}).*real(varargin{3}) + imag(varargin{1}).*imag(varargin{2}).*imag(varargin{3}) + ...
                        1i*(imag(varargin{1}).*real(varargin{2}).*real(varargin{3}) - real(varargin{1}).*imag(varargin{2}).*imag(varargin{3}));
                    self.z2 = real(varargin{1}).*imag(varargin{2}).*real(varargin{3}) - imag(varargin{1}).*real(varargin{2}).*imag(varargin{3}) + ...
                        1i*(real(varargin{1}).*real(varargin{2}).*imag(varargin{3}) + imag(varargin{1}).*imag(varargin{2}).*real(varargin{3}));
                otherwise
                    error('Requires exactly 2 or 3 inputs');
            end;
        end
        
%         function self = bicomplex(beta1,beta2)    % Idempotent representation
%             if nargin ~= 2
%                 error('Requires exactly 2 inputs')
%             end
%             if ~isequal(size(beta1),size(beta2))
%                 error('Inputs must be equally sized')
%             end
%             self.z1 = (beta1+beta2)/2;
%             self.z2 = 1i*(beta1-beta2)/2;
%         end
    end
    
    
    methods % Basic operators
        function mat = matrep(self) % Returns matrix representation
            mat = [self.z1,-self.z2;self.z2,self.z1];
        end
        
        function display(self)
            disp('z1:')
            disp(self.z1)
            disp('z2:')
            disp(self.z2)
        end
        
        function out = subsref(self,index) % Indexing
            if strcmp('()',index.type)
                out = bicomplex([],[]);
                out.z1 = builtin('subsref',self.z1,index);
                out.z2 = builtin('subsref',self.z2,index);
            elseif strcmp('.',index.type)
                out = eval(['self.',index.subs]);
            end
        end;
        
        function out = subsasgn_backup(self,index,value) % Asigning (BACKUP)
%             if strcmp('()',index.type)
%                 out = bicomplex([],[]);
%                 out.z1 = builtin('subsasgn',self.z1,index,value.z1);
%                 out.z2 = builtin('subsasgn',self.z2,index,value.z2);
%             elseif strcmp('.',index.type)
%                 if ~(strcmp(index.subs,'z1') || strcmp(index.subs,'z2'))
%                     error('No such field exists. Use z1 and z2 instead')
%                 else
%                     if strcmp(index.subs,'z1')
%                         z_1 = value;
%                         z_2 = self.z2;
%                     else
%                         z_2 = value;
%                         z_1 = self.z1;
%                     end
%                     out = bicomplex(z_1,z_2);
%                 end
%             end
        end
        
        function out = subsasgn(self,index,value) % Asigning
            if strcmp('()',index.type)
                z_1 = self.z1;
                z_2 = self.z2;
                z_1 = subsasgn(z_1, index, value.z1);            
                z_2 = subsasgn(z_2, index, value.z2);
            elseif strcmp('.',index.type)
                if ~(strcmp(index.subs,'z1') || strcmp(index.subs,'z2'))
                    error('No such field exists. Use z1 and z2 instead')
                    return;
                else
                    if strcmp(index.subs,'z1')
                        z_1 = value;
                        z_2 = self.z2;
                    else
                        z_2 = value;
                        z_1 = self.z1;
                    end
                end
            end
            out = bicomplex(z_1, z_2);
        end
        
        function out = horzcat(self,varargin) % Horizontal concatenation
            z_1 = [self.z1];
            z_2 = [self.z2];
            for i = 1:length(varargin)
                [~,tmp] = isbicomp([],varargin{i});
                z_1 = [z_1,tmp.z1];
                z_2 = [z_2,tmp.z2];
            end
            out = bicomplex(z_1,z_2);
        end
        
        function out = vertcat(self,varargin) % Vertical concatenation
            z_1 = [self.z1];
            z_2 = [self.z2];
            for i = 1:length(varargin)
                [~,tmp] = isbicomp([],varargin{i});
                z_1 = [z_1;tmp.z1];
                z_2 = [z_2;tmp.z2];
            end
            out = bicomplex(z_1,z_2);
        end
        
        function out = plus(self,other) % Addition
            [self,other] = isbicomp(self,other);
            if or(prod(size(self)), prod(size(other)))
                z_1 = self.z1 + other.z1;
                z_2 = self.z2 + other.z2;
                out = bicomplex(z_1, z_2);
            else
                zeta = matrep(self)+matrep(other);
                out = mat2bicomp(zeta);
            end;
        end
        
        function out = minus(self,other) % Subtraction
            [self,other] = isbicomp(self,other);
            zeta = matrep(self)- matrep(other);
            out = mat2bicomp(zeta);
        end
        
        function out = uplus(self) % Unary plus
            out = self;
        end
        
        function out = uminus(self) % Unary minus
            out = -1*self;
        end
        
%         function out = mtimes(self,other) % Multiplication
%             [self,other] = isbicomp(self,other);
%             if ~prod(size(self)==size(other)) && numel(self) == 1
%                 mat = matrep(self.*other);
%             elseif ~prod(size(self)==size(other)) && numel(other) == 1
%                 mat = matrep(self.*other);
%             else
%                 mat = matrep(self)*matrep(other);
%             end
%             out = mat2bicomp(mat);
%         end
        
        function out = mtimes(self,other) % Multiplication
            [self,other] = isbicomp(self,other);
            if ~all(size(self)==size(other)) && numel(self) == 1
                out = self.*other;
            elseif ~all(size(self)==size(other)) && numel(other) == 1
                out = self.*other;
            else
                z_1 = self.z1 * other.z1 - self.z2 * other.z2;
                z_2 = self.z1 * other.z2 + self.z2 * other.z1;
                out = bicomplex(z_1, z_2);
            end
        end

        function out = times(self,other) % Elementwise multiplication
            if and(numel(other)==1, isreal(other))
                z_1 = self.z1 * other;
                z_2 = self.z2 * other;
            elseif and(numel(self)==1, isreal(self))
                z_1 = other.z1 * self;
                z_2 = other.z2 * self;
            else
                [self,other] = isbicomp(self,other);
%                 if size(self) == size(other)
%                     sizes = size(self);
%                     z_1 = zeros(sizes);
%                     z_2 = zeros(sizes);
%                     for i = 1:prod(sizes)
%                         sr.type = '()';
%                         sr.subs = {i};
%                         tmp = subsref(self,sr)*subsref(other,sr);
%                         z_1(i) = tmp.z1;
%                         z_2(i) = tmp.z2;
%                     end

                if or(numel(self) == 1, numel(other) == 1)
                    z_1 = self.z1 .* other.z1 - self.z2 .* other.z2;
                    z_2 = self.z1 .* other.z2 + self.z2 .* other.z1;
                elseif size(self) == size(other)       % Need to cke it separately to allow multiplication of 3D arrays by a scalar
                    z_1 = self.z1 .* other.z1 - self.z2 .* other.z2;
                    z_2 = self.z1 .* other.z2 + self.z2 .* other.z1;
%                 elseif numel(self) == 1
%                     sizes = size(other);
%                     z_1 = zeros(sizes);
%                     z_2 = zeros(sizes);
%                     for i = 1:prod(sizes)
%                         sr.type = '()';
%                         sr.subs = {i};
%                         tmp = self*subsref(other,sr);
%                         z_1(i) = tmp.z1;
%                         z_2(i) = tmp.z2;
%                     end
%                 elseif numel(other) == 1
%                     sizes = size(self);
%                     z_1 = zeros(sizes);
%                     z_2 = zeros(sizes);
%                     for i = 1:prod(sizes)
%                         sr.type = '()';
%                         sr.subs = {i};
%                         tmp = subsref(self,sr)*other;
%                         z_1(i) = tmp.z1;
%                         z_2(i) = tmp.z2;
%                     end
                else
                    error('Matrix dimensions must agree')
                end
            end
            out = bicomplex(z_1,z_2);
        end
        
        function out = sum(self,varargin) % Summation
            if numel(varargin) == 1
                out = bicomplex(sum(self.z1,varargin{1}), sum(self.z2,varargin{1}));
            else
                out = bicomplex(sum(self.z1), sum(self.z2));
            end;
        end
        
        function out = mrdivide(self,other) % Division
            if numel(other) == 1 && numel(other) ~=numel(self)
                mat = matrep(self./other);
            else
                [self,other] = isbicomp(self,other);
                mat = matrep(self)/matrep(other);
            end
            out = mat2bicomp(mat);
        end
        
        function out = rdivide(self,other) % Elementwise division
            if and(numel(other)==1, isreal(other))
                z_1 = self.z1 / other;
                z_2 = self.z2 / other;
            else
                [self,other] = isbicomp(self,other);
                if size(self) == size(other)
                    sizes = size(self);
                    z_1 = zeros(sizes);
                    z_2 = zeros(sizes);
                    for i = 1:prod(sizes)
                        
                        sr.type = '()';
                        sr.subs = {i};
                        tmp = subsref(self,sr)/subsref(other,sr);
                        z_1(i) = tmp.z1;
                        z_2(i) = tmp.z2;
                    end
                elseif numel(self) == 1
                    sizes = size(other);
                    z_1 = zeros(sizes);
                    z_2 = zeros(sizes);
                    for i = 1:prod(sizes)
                        sr.type = '()';
                        sr.subs = {i};
                        tmp = self/subsref(other,sr);
                        z_1(i) = tmp.z1;
                        z_2(i) = tmp.z2;
                    end
                elseif numel(other) == 1
                    sizes = size(self);
                    z_1 = zeros(sizes);
                    z_2 = zeros(sizes);
                    for i = 1:prod(sizes)
                        sr.type = '()';
                        sr.subs = {i};
                        tmp = subsref(self,sr)/other;
                        z_1(i) = tmp.z1;
                        z_2(i) = tmp.z2;
                    end
                else
                    error('Matrix dimensions must agree.')
                end
            end;
            out = bicomplex(z_1,z_2);
        end
        
        function out = inv(self) % Matrix inverse
            sizes = size(self);
            if size(1) == size(2)
                out = eye(sizes(1)) / self;
            else
                error('Matrix should be square');
            end;
        end
        
        function out = power(self,other) % Elementwise power
            sizes = size(self);
            z_1 = zeros(sizes);
            z_2 = zeros(sizes);
            
            for i = 1:length(z_1(:))
                sr.type = '()';
                sr.subs = {i};
                r = modc(subsref(self,sr));
                theta = argc(subsref(self,sr));
                z_1(i) = r^other*cos(other*theta);
                z_2(i) = r^other*sin(other*theta);
            end
            out = bicomplex([],[]);
            out.z1 = z_1;
            out.z2 = z_2;
            
        end
        
        function out = mpower(self,other) % Elementwise power
            sizes = size(self);
            z_1 = zeros(sizes);
            z_2 = zeros(sizes);
            
            for i = 1:length(z_1(:))
                sr.type = '()';
                sr.subs = {i};
                r = modc(subsref(self,sr));
                theta = argc(subsref(self,sr));
                z_1 = r^other*cos(other*theta);
                z_2 = r^other*sin(other*theta);
            end
            out = bicomplex([],[]);
            out.z1 = z_1;
            out.z2 = z_2;
            
        end
        
        function dims = size(self) % Returning size of array
            dims = size(self.z1);
        end
        
        function n = numel(self) % Returning number of elements
            n = numel(self.z1);
        end
        
        function out = modc(self) % Complex modulus
            out = sqrt(self.z1.^2 + self.z2.^2);
        end
        
        function out = norm(self) % Norm
            out = sqrt(real(self.z1).^2 + real(self.z2).^2 + ...
                imag(self.z1).^2 + imag(self.z2).^2);
        end
        
        function out = det(self) % Matrix determinant
            sizes = size(self);
            if sizes(1) == sizes(2)
                out = det(self.z1 - 1i*self.z2) + det(self.z1 + 1i*self.z2);
            else
                error('Matrix should be square');
            end;
        end
        
        function theta = argc(self) % Complex argument
            theta = atan2(self);
        end
        
        function out = conj(self) % Bicomplex conjugate (both i and j are conjugated)
            out = self;
            out.z1 = conj(out.z1);
            out.z2 = - conj(out.z2);
        end
        
        function out = conji(self) % Ci-conjugate
            out = self;
            out.z1 = conj(out.z1);
            out.z2 = conj(out.z2);
        end
        
        function out = conjj(self) % Cj-conjugate
            out = self;
            out.z2 = - out.z2;
        end
        
        function out = transpose(self) % Transposition (without the conjugate)
            z_1 = transpose(self.z1);
            z_2 = transpose(self.z2);
            out = bicomplex(z_1, z_2);
        end
        
        function out = diag(self)   % Diagonal of a matrix
            sz = size(self);
            if sz(1) == sz(2)
                out = bicomplex(diag(real(self)) + 1i*diag(imag1(self)), diag(imag2(self)) + 1i*diag(imag12(self)) );
            end
        end
        
        function out = reshape(self, varargin) % Reshape
            if nargin == 2+1
                z_1 = reshape(self.z1, varargin{1}, varargin{2});
                z_2 = reshape(self.z2, varargin{1}, varargin{2});
            elseif nargin == 3+1
                z_1 = reshape(self.z1, varargin{1}, varargin{2}, varargin{3});
                z_2 = reshape(self.z2, varargin{1}, varargin{2}, varargin{3});
            end;
            out = bicomplex(z_1, z_2);
        end
        
        function out = repmat(self, varargin) % Reshape
            if nargin == 2+1
                z_1 = repmat(self.z1, varargin{1}, varargin{2});
                z_2 = repmat(self.z2, varargin{1}, varargin{2});
            elseif nargin == 3+1
                z_1 = repmat(self.z1, varargin{1}, varargin{2}, varargin{3});
                z_2 = repmat(self.z2, varargin{1}, varargin{2}, varargin{3});
            end;
            out = bicomplex(z_1, z_2);
        end
        
        function out = lt(self,other) % Less than, self < other
            out = false;
            if real(self.z1) < real(other.z1)
                out = true;
            end
        end
        
        function out = gt(self,other) % Greater than, self > other
            out = false;
            if real(self.z1) > real(other.z1)
                out = true;
            end
        end
        
        function out = le(self,other) % Less than or equal, self <= other
            out = false;
            if real(self.z1) <= real(other.z1)
                out = true;
            end
        end
        
        function out = ge(self,other) % Greater than or equal, self >= other
            out = false;
            if real(self.z1) >= real(other.z1)
                out = true;
            end
        end
        
        function out = eq(self,other) % Equality, self == other
            out = false;
            if self.z1 == other.z1 && self.z2 == other.z2
                out = true;
            end
        end
        
        function out = ne(self,other) % Not equal, self ~= other
            out = true;
            if self.z1 == other.z1 && self.z2 == other.z2
                out = false;
            end
        end
    end
    
    methods % Mathematical functions
        %% Exponential function and logarithm
        function out = exp(self) % Exponential
            
            
            out = bicomplex([],[]);
            out.z1=exp(self.z1).*cos(self.z2);
            out.z2=exp(self.z1).*sin(self.z2);
        end
        
        function out = log(self) % Natural logaritm
            out = bicomplex([],[]);
            out.z1=log(modc(self));
            out.z2=argc(self);
        end
        
        %% Basic trigonometric functions
        function out = sin(self) % sin
            out = bicomplex([],[]);
            out.z1=cosh(self.z2).*sin(self.z1);
            out.z2=sinh(self.z2).*cos(self.z1);
        end
        
        function out = cos(self) % cos
            out = bicomplex([],[]);
            out.z1=cosh(self.z2).*cos(self.z1);
            out.z2=-sinh(self.z2).*sin(self.z1);
        end
        
        function out = tan(self) % tan
            out = sin(self)./cos(self);
        end
        
        function out = cot(self) % cot
            out = cos(self)./sin(self);
        end
        
        function out = sec(self) % sec
            out = 1./cos(self);
        end
        
        function out = csc(self) % csc
            out = 1./sin(self);
        end
        
        %% Basic hyperbolic functions
        function out = sinh(self)
            out = bicomplex([],[]);
            out.z1=cosh(self.z1).*cos(self.z2);
            out.z2=sinh(self.z1).*sin(self.z2);
        end
        
        function out = cosh(self)
            out = bicomplex([],[]);
            out.z1=sinh(self.z1).*cos(self.z2);
            out.z2=cosh(self.z1).*sin(self.z2);
        end
        
        function out = tanh(self)
            
            out = sinh(self)./cosh(self);
        end
        
        function out = coth(self)
            out = cosh(self)./sinh(self);
        end
        
        function out = sech(self)
            out = 1./cosh(self);
        end
        
        function out = csch(self)
            out = 1./sinh(self);
        end
        
        function out = atan2(self)
            sizes = size(self);
            ang = zeros(sizes);
            
            for i = 1:prod(sizes)
                sr.type = '()';
                sr.subs = {i};
                if real(self.z1(i)) > 0;
                    ang(i) = atan(self.z2(i)./ self.z1(i));
                elseif real(self.z1(i))<0 && real(self.z2(i))>= 0;
                    ang(i) = atan(self.z2(i)./self.z1(i))+pi;
                elseif real(self.z1(i))<0 && real(self.z2(i))<0;
                    ang(i) = atan(self.z2(i)./self.z1(i))-pi;
                elseif real(self.z1(i))==0 && real(self.z2(i))> 0;
                    ang(i) = pi/2;
                elseif real(self.z1(i))==0 && real(self.z2(i))< 0;
                    ang(i) = -pi/2;
                else
                    error('atan(0,0) undefined');
                end
            end
            out = ang;
        end
        
        function out = sqrt(self)
            out = self.^0.5;
        end
        
%         %% Fourier transform
%         function out = fft(self, varargin) % Exponential
%             nt = size(self); nt1 = nt(1); nt2 = nt(2);
%             out = bicomplex([],[]);
%             out.z1=exp(self.z1).*cos(self.z2);
%             out.z2=exp(self.z1).*sin(self.z2);
%         end
    end
    
    methods % Functions for returning the imaginary and real parts
        function out = cmpl1(self)
            out = self.z1;
        end;
        function out = cmpl2(self)
            out = self.z2;
        end;
        function out = real(self)
            out = real(self.z1);
        end
        function out = imag1(self)
            out = imag(self.z1);
        end
        function out = imag2(self)
            out = real(self.z2);
        end
        function out = imag12(self)
            out = imag(self.z2);
        end
        function [out1, out2, out3] = angles(self)
            % Returns the three argumanets in a polar representation, exp(1i*out1)*exp(1j*out2)*exp(1k*out3)
            d = sqrt(self .* conj(self));
            self = self ./ d;
            out1 = angle(self.z1);
            out2 = angle(real(self.z1) + 1i*real(self.z2));
            out3 = angle(real(d.z1) + 1i*imag(d.z2));
        end
    end
    
    methods %% Utility functions
        function [self,other] = isbicomp(self,other)
            % Verifies that self and other are bicomplex, or converts them to bicomplex if possible
            if isa(self,'double')
                self = bicomplex(self,zeros(size(self)));
            elseif ~isa(self,'bicomplex')
                error('Self is not of class bicomplex')
            end
            if isa(other,'double')
                other = bicomplex(other,zeros(size(other)));
            elseif ~isa(other,'bicomplex')
                error('Other is not of class bicomplex')
            end
        end
        
        function [out] = snr(self, other)
            out = snr([cmpl1(self); cmpl2(self)], [cmpl1(other); cmpl2(other)]);
        end
        
%         function zeta = mat2bicomp(mat)
%             % Takes the matrix representation and returns the corresponding bicomplex
%             sizes = size(mat);
%             str1 = '1:sizes(1)/2,1:sizes(2)/2';
%             str2 = 'sizes(1)/2+1:end,1:sizes(2)/2';
%             for i = 3:length(sizes);
%                 str1 = [str1 sprintf(',1:sizes(%i)',i)];
%                 str2 = [str2 sprintf(',1:sizes(%i)',i)];
%             end
%             str1 = sprintf('mat(%s)',str1);
%             str2 = sprintf('mat(%s)',str2);
%             z1 = eval(str1);
%             z2 = eval(str2);
%             zeta = bicomplex(z1, z2);
%         end
    end
    
end

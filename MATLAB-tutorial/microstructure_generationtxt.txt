clear
clc
close all

el=101;  % Size of the microstructure

nods=50; % Number of data sets per class

% Isotropic datasets
M=zeros(nods,el^2);
[X,Y]=meshgrid(-floor(el/2):floor(el/2));
F=mvnpdf([X(:) Y(:)],[0 0],[10 0; 0 10]);
F2=reshape(F,el,el);

for ii=1:nods    
    M2=ifft2(fft2(F2).*conj(fft2(rand(el,el))));
    M2s=sort(M2(:));
    M2(M2>M2s(round(ceil(el^2/2)+randn(1)*300)))=1;
    M2(M2~=1)=0;
    M(ii,:)=M2(:);    
end

save M1.mat M

% X fiber datasets
M=zeros(nods,el^2);
[X,Y]=meshgrid(-floor(el/2):floor(el/2));
F=mvnpdf([X(:) Y(:)],[0 0],[20 0; 0 2]);
F2=reshape(F,el,el);

for ii=1:nods    
    M2=ifft2(fft2(F2).*conj(fft2(rand(el,el))));
    M2s=sort(M2(:));
    M2(M2>M2s(round(ceil(el^2/2)+randn(1)*300)))=1;
    M2(M2~=1)=0;
    M(ii,:)=M2(:);    
end

save M2.mat M

% Y fiber datasets
M=zeros(nods,el^2);
[X,Y]=meshgrid(-floor(el/2):floor(el/2));
F=mvnpdf([X(:) Y(:)],[0 0],[2 0; 0 20]);
F2=reshape(F,el,el);

for ii=1:nods    
    M2=ifft2(fft2(F2).*conj(fft2(rand(el,el))));
    M2s=sort(M2(:));
    M2(M2>M2s(round(ceil(el^2/2)+randn(1)*300)))=1;
    M2(M2~=1)=0;
    M(ii,:)=M2(:);    
end

save M3.mat M
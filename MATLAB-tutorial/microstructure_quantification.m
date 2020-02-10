%% Microstructure Classification tutorial
%
% This script contains the fundamental codes for calculation of two-point
% statistics and principal component analysis to obtain reduced order
% statistical representation of microstructure.
%
% All necessary files along with a set of slides with additional
% information will be included in MATIN platform.
%
% For questions about the tutorial, you can send me an email.
% yabansu@gatech.edu
%
% If you are using (or planing to use) MATLAB for your class project and if
% you need help with using it, you can schedule a meeting with me to
% discuss how you can use MATLAB effectively for materials informatics
% applications.
%
% September 25, 2017
% Yuksel C. Yabansu, PhD

clear
clc
close all

% Preparing a full raw data matrix for the ensemble of 3 microstructure
% classes
Mf=[];  % Initialization of data matrix
load M1.mat
Mf=[Mf; M];
load M2.mat
Mf=[Mf; M];
load M3.mat
Mf=[Mf; M];

% Variables
noms=size(Mf,1)/3;  % Number of microstructures in each class

el=round(sqrt(size(M,2)));  % Microstructure size

nods=size(M,1);  % Number of datasets in each class of microstructures

% Truncation size for two-point statistics. This number includes the center
% cell which is r=0 vector. For example, trunc=31 cuts off the center
% region with a size of 61x61 from the full statistics map. This box is
% specified by center cell (1) plus maximum vector size of 30 (in total it
% makes 31).
trunc=21;

%% Visualization of microstructures
% Each row of the microstructure matrix is one microstructure (or one
% observation). For the purpose of visualization of microstructure and its
% two-point statistics, we have to convert it to 2-D matrix. We use reshape
% function to convert a row vector to array.
M1=reshape(Mf(1,:),el,el);  % Microstructure with isotropic morphology
M2=reshape(Mf(nods+1,:),el,el);  % Microstructure with isotropic morphology
M3=reshape(Mf(2*nods+1,:),el,el);  % Microstructure with isotropic morphology

% Visualization of sample microstructure from each class
figure
subplot(1,3,1)
imagesc(M1) % visualize the Image
axis image % preserves the aspect ratio
colormap('gray') % specifies the grey larges value white other 0 , if we woul make it jet it would make it different!!!
set(gca,'FontSize',16)
title({'Class 1','Isotropic morphology'},'FontSize',16)

subplot(1,3,2)
imagesc(M2)
axis image
colormap('gray')
set(gca,'FontSize',16)
title({'Class 2','Oriented in x'},'FontSize',16)

subplot(1,3,3)
imagesc(M3)
axis image
colormap('gray')
set(gca,'FontSize',16)
title({'Class 3','Oriented in y'},'FontSize',16)

%% Calculation of two-point statistics
% For a material system with H distinct local states, there are H^2 sets of
% two-point statistics. For the two-phase materials system, there 4 sets,
% which 2 are autocorrelations and other 2 are cross-correlation.

% Two-point autocorrelation is defined by setting the same local state at
% the head and tail of the vectors (i.e. h=h'). Crosscorrelation is when
% there are different local states at the head and tail of the vector (i.e.
% h~=h')

% Autocorrelation can be calculated as a convolution of microstructure
% signal with itself. In discrete Fourier space, this convolution is
% nothing but a element-by-element multiplication. Crosscorrelation is
% calculated the same way except the microstructure signal is convolved
% with the signal of other phase.

% numel(M1) is the normalization factor. For periodic microstructure, it is
% a constant for every possible vector. For nonperiodic microstructures, a
% separate normalization matrix has to be arranged.

% Full statistics of microstructures
F11=fftshift(ifft2(fft2(M1).*conj(fft2(M1))))/numel(M1);
F00=fftshift(ifft2(fft2(M1==0).*conj(fft2(M1==0))))/numel(M1);
F10=fftshift(ifft2(fft2(M1).*conj(fft2(M1==0))))/numel(M1);
F01=fftshift(ifft2(fft2(M1==0).*conj(fft2(M1))))/numel(M1);

% Visualization of two-point statistics of example microstructure from each
% class
% Please pay attention to the axes when plotting two-point statistics.
% Vector space that two-point statistics is defined has the same grid as
% the microstructure but the coordinate system is different.

figure
ax1=subplot(2,4,[1 2 5 6]);
imagesc(M1)
axis image
colormap(ax1,'gray')
set(gca,'FontSize',16)
title('Microstructure','FontSize',16)

subplot(2,4,3)
imagesc(-floor(el/2),-floor(el/2),F00)
axis image
colormap('jet')
set(gca,'FontSize',16)
xlabel('r_x','Fontsize',16,'FontWeight','bold')
ylabel('r_y','Fontsize',16,'FontWeight','bold')
title('f_{r}^{00}','FontSize',16)
cb=colorbar;
set(cb,'FontSize',16)

subplot(2,4,4)
imagesc(-floor(el/2),-floor(el/2),F01)
axis image
colormap('jet')
set(gca,'FontSize',16)
xlabel('r_x','Fontsize',16,'FontWeight','bold')
ylabel('r_y','Fontsize',16,'FontWeight','bold')
title('f_{r}^{01}','FontSize',16)
cb=colorbar;
set(cb,'FontSize',16)

subplot(2,4,7)
imagesc(-floor(el/2),-floor(el/2),F10)
axis image
colormap('jet')
set(gca,'FontSize',16)
xlabel('r_x','Fontsize',16,'FontWeight','bold')
ylabel('r_y','Fontsize',16,'FontWeight','bold')
title('f_{r}^{10}','FontSize',16)
cb=colorbar;
set(cb,'FontSize',16)

subplot(2,4,8)
imagesc(-floor(el/2),-floor(el/2),F11)
axis image
colormap('jet')
set(gca,'FontSize',16)
xlabel('r_x','Fontsize',16,'FontWeight','bold')
ylabel('r_y','Fontsize',16,'FontWeight','bold')
title('f_{r}^{11}','FontSize',16)
cb=colorbar;
set(cb,'FontSize',16)

%% Visualization of two-point statistics for all 3 classes
% Full statistics of microstructures
F1=fftshift(ifft2(fft2(M1).*conj(fft2(M1))))/numel(M1);
F2=fftshift(ifft2(fft2(M2).*conj(fft2(M2))))/numel(M2);
F3=fftshift(ifft2(fft2(M3).*conj(fft2(M3))))/numel(M3);

% Truncating the statistics
F1t=F1(ceil(el/2)-trunc+1:ceil(el/2)+trunc-1,...
    ceil(el/2)-trunc+1:ceil(el/2)+trunc-1);
F2t=F2(ceil(el/2)-trunc+1:ceil(el/2)+trunc-1,...
    ceil(el/2)-trunc+1:ceil(el/2)+trunc-1);
F3t=F3(ceil(el/2)-trunc+1:ceil(el/2)+trunc-1,...
    ceil(el/2)-trunc+1:ceil(el/2)+trunc-1);

% Visualization of two-point statistics of example microstructure from each
% class
% Please pay attention to the axes when plotting two-point statistics.
% Vector space that two-point statistics is defined has the same grid as
% the microstructure but the coordinate system is different.
figure
subplot(2,3,1)
imagesc(-floor(el/2),-floor(el/2),F1)
axis image
colormap('jet')
set(gca,'FontSize',16)
xlabel('r_x','Fontsize',16,'FontWeight','bold')
ylabel('r_y','Fontsize',16,'FontWeight','bold')
title({'f_{r}^{00}'},'FontSize',16)
cb=colorbar;
set(cb,'FontSize',16)

subplot(2,3,2)
imagesc(-floor(el/2),-floor(el/2),F2)
axis image
colormap('jet')
set(gca,'FontSize',16)
xlabel('r_x','Fontsize',16,'FontWeight','bold')
ylabel('r_y','Fontsize',16,'FontWeight','bold')
title({'Class 2','Full f_{r}^{11}'},'FontSize',16)
cb=colorbar;
set(cb,'FontSize',16)

subplot(2,3,3)
imagesc(-floor(el/2),-floor(el/2),F3)
axis image
colormap('jet')
set(gca,'FontSize',16)
xlabel('r_x','Fontsize',16,'FontWeight','bold')
ylabel('r_y','Fontsize',16,'FontWeight','bold')
title({'Class 3','Full f_{r}^{11}'},'FontSize',16)
cb=colorbar;
set(cb,'FontSize',16)

subplot(2,3,4)
imagesc(-trunc+1,-trunc+1,F1t)
axis image
colormap('jet')
set(gca,'FontSize',16)
xlabel('r_x','Fontsize',16,'FontWeight','bold')
ylabel('r_y','Fontsize',16,'FontWeight','bold')
title('Truncated f_{r}^{11}','FontSize',16)



subplot(2,3,5)
imagesc(-trunc+1,-trunc+1,F2t)
axis image
colormap('jet')
set(gca,'FontSize',16)
xlabel('r_x','Fontsize',16,'FontWeight','bold')
ylabel('r_y','Fontsize',16,'FontWeight','bold')
title('Truncated f_{r}^{11}','FontSize',16)

subplot(2,3,6)
imagesc(-trunc+1,-trunc+1,F3t)
axis image
colormap('jet')
set(gca,'FontSize',16)
xlabel('r_x','Fontsize',16,'FontWeight','bold')
ylabel('r_y','Fontsize',16,'FontWeight','bold')
title('Truncated f_{r}^{11}','FontSize',16)

%% Principal component analysis for microstructure classification
% Initialization of data matrix. This matrix will go in PCA function. The
% rows are observations which are two-point autocorrelations of white
% phases in microstructures. The columns are the number of statistics in
% the truncated autocorrelations (i.e. (2*trunc-1)^2)
F=zeros(size(Mf,1),(2*trunc-1)^2);
for ii=1:size(Mf,1)
    M2=reshape(Mf(ii,:),el,el);
    F1=fftshift(ifft2(fft2(M2).*conj(fft2(M2))))/numel(M2);
    F1t=F1(ceil(el/2)-trunc+1:ceil(el/2)+trunc-1,...
        ceil(el/2)-trunc+1:ceil(el/2)+trunc-1);
    F(ii,:)=F1t(:);
end

% Mean of the data
F_mean=mean(F,1);

% Principal component analysis
[coeff,score,~,~,explained]=pca(F,'NumComponents',3);

figure
cmap=[1 0 0; 0 1 0; 0 0 1];
subplot(2,2,1)
for ii=1:3
    plot3(score((ii-1)*noms+1:ii*noms,1),score((ii-1)*noms+1:ii*noms,2),...
        score((ii-1)*noms+1:ii*noms,3),'ko','MarkerFaceColor',cmap(ii,:),...
        'MarkerSize',8)
    hold on
end
grid on
axis image
set(gca,'FontSize',16)
xlabel('PC1','FontSize',16)
ylabel('PC2','FontSize',16)
zlabel('PC3','FontSize',16)
l=legend('Class 1 - Isotropic','Class 2 - X Oriented','Class 3 - Y Oriented');
set(l,'FontSize',16)

subplot(2,2,2)
for ii=1:3
    plot(score((ii-1)*noms+1:ii*noms,1),score((ii-1)*noms+1:ii*noms,2),...
        'ko','MarkerFaceColor',cmap(ii,:),'MarkerSize',8)
    hold on
end
grid on
set(gca,'FontSize',16)
xlabel('PC1','FontSize',16)
ylabel('PC2','FontSize',16)
xlim([min(score(:,1)) max(score(:,1))])
ylim([min(score(:,2)) max(score(:,2))])

subplot(2,2,3)
for ii=1:3
    plot(score((ii-1)*noms+1:ii*noms,1),score((ii-1)*noms+1:ii*noms,3),...
        'ko','MarkerFaceColor',cmap(ii,:),'MarkerSize',8)
    hold on
end
grid on
set(gca,'FontSize',16)
xlabel('PC1','FontSize',16)
ylabel('PC3','FontSize',16)
xlim([min(score(:,1)) max(score(:,1))])
ylim([min(score(:,3)) max(score(:,3))])

subplot(2,2,4)
for ii=1:3
    plot(score((ii-1)*noms+1:ii*noms,2),score((ii-1)*noms+1:ii*noms,3),...
        'ko','MarkerFaceColor',cmap(ii,:),'MarkerSize',8)
    hold on
end
grid on
set(gca,'FontSize',16)
xlabel('PC2','FontSize',16)
ylabel('PC3','FontSize',16)
xlim([min(score(:,2)) max(score(:,2))])
ylim([min(score(:,3)) max(score(:,3))])

%% Visualization of explained variance

figure
pareto(explained)
grid on
xlabel('Principal Component','FontSize',16)
ylabel('Explained variance (%)','FontSize',16)
set(gca,'FontSize',16)

%% Visualization of basis vectors
% Basis vectors are stored as column vectors in variable coeff.
phi1=reshape(coeff(:,1),2*trunc-1,2*trunc-1);
phi2=reshape(coeff(:,2),2*trunc-1,2*trunc-1);
phi3=reshape(coeff(:,3),2*trunc-1,2*trunc-1);

% 2-D matrix of average statistics
F_mean_2D=reshape(F_mean,2*trunc-1,2*trunc-1);

figure
subplot(2,2,1)
imagesc(-trunc+1,-trunc+1,F_mean_2D)
axis image
colormap('jet')
xlabel('r_x','Fontsize',16,'FontWeight','bold')
ylabel('r_y','Fontsize',16,'FontWeight','bold')
title('Average statistics','FontSize',16)
cb=colorbar;
set(cb,'FontSize',16)
set(gca,'FontSize',16)

subplot(2,2,2)
imagesc(-trunc+1,-trunc+1,phi1)
axis image
colormap('jet')
xlabel('r_x','Fontsize',16,'FontWeight','bold')
ylabel('r_y','Fontsize',16,'FontWeight','bold')
title('\phi_{1,r}^{11}','FontSize',16)
cb=colorbar;
set(cb,'FontSize',16)
set(gca,'FontSize',16)
caxis([min(coeff(:)) max(coeff(:))])

subplot(2,2,3)
imagesc(-trunc+1,-trunc+1,phi2)
axis image
colormap('jet')
xlabel('r_x','Fontsize',16,'FontWeight','bold')
ylabel('r_y','Fontsize',16,'FontWeight','bold')
title('\phi_{2,r}^{11}','FontSize',16)
cb=colorbar;
set(cb,'FontSize',16)
set(gca,'FontSize',16)
caxis([min(coeff(:)) max(coeff(:))])

subplot(2,2,4)
imagesc(-trunc+1,-trunc+1,phi3)
axis image
colormap('jet')
xlabel('r_x','Fontsize',16,'FontWeight','bold')
ylabel('r_y','Fontsize',16,'FontWeight','bold')
title('\phi_{3,r}^{11}','FontSize',16)
cb=colorbar;
set(cb,'FontSize',16)
set(gca,'FontSize',16)
caxis([min(coeff(:)) max(coeff(:))])
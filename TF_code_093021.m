% Thin film transmission code for 4 layers
% Allows user to plot theoretical transmission and phase plots for
% different thicknesses of thin-films
%Sanket Deshpande

clear all;
close all;
digits(32)

%% Import the data and convert the wavenumber to wavelength (in ascending order)

dataImport = readmatrix('G:\My Drive\Projects\Quantum Array Generator\FTIR Data\T_072721_TiO2_310_SiO2.dpt', 'FileType', 'text');
data = zeros(size(dataImport));

for i=1:size(dataImport,1)
    data(i,1)=1e7/dataImport(size(dataImport,1)-(i-1),1);
    data(i,2)=dataImport(size(dataImport,1)-(i-1),2);
end

%% Set up the simulation

%wavelength range. Use this if you're defining the refractive indices using Sellmeier relation
%lambda = linspace(data(1,1),data(size(data,1),1),size(data,1)/10);

%use this definition of wavelength if importing refractive indices of
%materials from CSV file

TiO2_ref = readmatrix('G:\My Drive\Projects\Quantum Array Generator\Thin Film matlab codes\TF theory vs data code\tio2_Ellipsometry_nk_300_1000.txt', 'FileType', 'text');
lambda(1,:) = TiO2_ref(:,1);
% TiO2_ref = readmatrix('G:\My Drive\Projects\Quantum Array Generator\Thin Film matlab codes\TF theory vs data code\TiO2_n.csv', 'FileType', 'text');
% lambda(1,:) = TiO2_ref(:,1).*1000;
%SiO2_ref = readmatrix('G:\My Drive\Projects\Quantum Array Generator\Thin Film matlab codes\TF theory vs data code\SiO2_n.csv', 'FileType', 'text');


%% free-space refractive indices of various media

%n_FSi = 1.4531; %fused silica refractive index at 810 nm
% n_FSi = sqrt(1 + 0.6961663.*((lambda./1000).^2)./((lambda./1000).^2 - 0.0684043^2) + 0.4079426.*((lambda./1000).^2)./((lambda./1000).^2 - 0.1162414^2) + 0.8974794.*((lambda./1000).^2)./((lambda./1000).^2 - 9.896161^2) ); %Sellmeier relation from https://refractiveindex.info/?shelf=glass&book=fused_silica&page=Malitson
n_FSi = zeros(length(lambda),1);

% for j=1:size(SiO2_ref,1)
%     n_FSi(j,1) = SiO2_ref(j,2);
% end

%Cauchy parameters for FSi determined from ellipsometry measurements
A=1.4406;
B=0.0017069;
C=0.00012464;

for j=1:length(lambda)
    n_FSi(j,1) = A + B/((lambda(1,j)*0.001)^2) + C/((lambda(1,j)*0.001)^4); %Cauchy's equation. Lambda should be in micrometers
end    

n_air = ones(length(lambda),1); %air refractive index

%n_TiO2 = 2.5173; %Titanium dioxide at 810 nm

n_TiO2 = zeros(length(lambda),1);
for j=1:size(TiO2_ref,1)
    n_TiO2(j,1) = TiO2_ref(j,2);
end

%% variable declarations for different types of structures

nA = vpa(zeros(6,length(lambda))); %matrix holding refractive indices of all layers of Stucture A
lA = zeros(6,1); %matrix holding lengths of all layers of Structure A. Length units are nanometers

nB = vpa(zeros(6,length(lambda))); %matrix holding refractive indices of all layers of Stucture B
lB = zeros(6,1); %matrix holding lengths of all layers of Structure B. Length units are nanometers

nC = vpa(zeros(6,length(lambda))); %matrix holding refractive indices of all layers of Stucture C
lC = zeros(6,1); %matrix holding lengths of all layers of Structure C. Length units are nanometers

%% defining layers in the structure A
%if you intend to use less than 4 layers, you can define layers 2 and 3 to
%be made up of the same material. For example, 20 nm of Au can be defined
%as 10 nm Au in layer 1 and 10 nm Au in layer 2

%layer 1 is air with no length
nA(1,:) = n_air;
lA(1) = 0;

%layer 2 (thin film)
nA(2,:)=n_TiO2;
lA(2)=95;

%layer 3 (thin film)
nA(3,:)=n_FSi;
lA(3)=50;

%layer 4 (thin film)
nA(4,:)=n_FSi;
lA(4)=10;

%layer 5 (thin film)
nA(5,:)=n_FSi;
lA(5)=10;

%layer 6 (substrate)
nA(6,:)=n_FSi;
lA(6)=10;


%% defining layers in the structure B
%if you intend to use less than 4 layers, you can define layers 2 and 3 to
%be made up of the same material. For example, 20 nm of Au can be defined
%as 10 nm Au in layer 1 and 10 nm Au in layer 2

%layer 1 is air with no length
nB(1,:) = n_air;
lB(1) = 0;

%layer 2 (thin film)
nB(2,:)=n_TiO2;
lB(2)=100;

%layer 3 (thin film)
nB(3,:)=n_FSi;
lB(3)=50;

%layer 4 (thin film)
nB(4,:)=n_FSi;
lB(4)=5;

%layer 5 (thin film)
nB(5,:)=n_FSi;
lB(5)=10;

%layer 6 (substrate)
nB(6,:)=n_FSi;
lB(6)=10;

%% defining layers in the structure C
%if you intend to use less than 4 layers, you can define layers 2 and 3 to
%be made up of the same material. For example, 20 nm of Au can be defined
%as 10 nm Au in layer 1 and 10 nm Au in layer 2

%layer 1 is air with no length
nC(1,:) = n_air;
lC(1) = 0;

%layer 2 (thin film)
nC(2,:)=n_TiO2;
lC(2)=105;

%layer 3 (thin film)
nC(3,:)=n_FSi;
lC(3)=140;

%layer 4 (thin film)
nC(4,:)=n_FSi;
lC(4)=140;

%layer 5 (thin film)
nC(5,:)=n_FSi;
lC(5)=0;

%layer 6 (substrate)
nC(6,:)=n_FSi;
lC(6)=42;

%% Determining and plotting the transmissions and phase shifts for all structures

[TA,phaseA] = TS(nA,lA,lambda);
[TB,phaseB] = TS(nB,lB,lambda);
[TC,phaseC] = TS(nC,lC,lambda);
%plot the transmission vs lambda curve for TiO2 on SiO2
figure(1)
plot(lambda,TA.*100,'-','Color','g','Linewidth',3)
hold on
plot(lambda,TB.*100,'-','Color','r','Linewidth',3)
plot(lambda,TC.*100,'-','Color','b','Linewidth',3)
plot(data(:,1),data(:,2).*100,'*','Color','g','Linewidth',3)
ylabel('Transmittance [%]')
yyaxis right
%plot(lambda,phase,'-','Color','r','Linewidth',3)
xline(810,'-',{'810 nm'},'LabelVerticalAlignment','bottom','LabelOrientation','horizontal','Linewidth',2,'FontSize',25);
set(gca, 'XColor','k', 'YColor','k')
ylabel('Phase shift [deg]', 'Color','r')
hold off
title("TiO2 on SiO2")
xlabel('Free Space Wavelength [nm]')

%legend('T (theoretical)','T (experimental, IT=250)','\phi (theoretical)')
legend("t_{TiO2}="+lA(2,1),"t_{TiO2}="+lB(2,1),"t_{TiO2}="+lC(2,1))
ax = gca;
ax.FontSize = 45;

%% function to determine the transmission and phase shift for input parameters: refractive index sequence, thickness sequence, wavelength range
function [T,phase] = TS(n_mat,l_mat,lambda)
    t_mat = zeros(5,length(lambda)); %transmission coeff matrix for all interfaces
    r_mat = zeros(5,length(lambda)); %reflection coeff matrix for all interfaces
    phi_mat = zeros(6,length(lambda)); %phase shift matrix for all layers
    T_mat = zeros(2,2,5,length(lambda)); %Transfer Matrix for interfaces
    P_mat = zeros(2,2,6,length(lambda)); %Phase shift Matrix for layers
    M = zeros(2,2,length(lambda)); %Total transfer matrix
    for i=1:length(lambda)
        M(:,:,i) = [1 0;0 1];
    end
    
    %calculating the transfer matrix
    for i=1:length(lambda)
        for a=1:size(t_mat,1)
            t_mat(a,i)= 2*n_mat(a,i)/(n_mat(a,i)+n_mat(a+1,i));
            r_mat(a,i)= (n_mat(a,i)-n_mat(a+1,i))/(n_mat(a,i)+n_mat(a+1,i));

            T_mat(:,:,a,i) = t_mat(a,i)^(-1) * [1 r_mat(a,i); r_mat(a,i) 1];
        end
    end

    %calculating the propagation matrix
    for i=1:length(lambda)
        for b=1:(size(phi_mat,1)-1)
            phi_mat(b,i) = (2*pi/lambda(i))*n_mat(b+1,i)*l_mat(b+1);
            P_mat(:,:,b,i) = [exp(-1i*phi_mat(b,i)) 0; 0 exp(1i*phi_mat(b,i))];
        end
    end

    %calculating the total transfer matrix
    for i=1:length(lambda)
        for c=1:size(T_mat,1)
            M(:,:,i) = M(:,:,i) * T_mat(:,:,c,i) * P_mat(:,:,c,i);
        end
    end

    phase=zeros(length(lambda),1);
    T=zeros(length(lambda),1);

    for i=1:length(lambda)
        %phase shift due to all layers in degrees
        phase(i) = angle(1/M(1,1,i)).*180/pi;

        %Transmittance from all the layers
        T(i)=(n_mat(size(n_mat,1),i)/n_mat(1,i))*(abs(1/M(1,1,i)))^2;
    end
end




%% Code reservoir

% %plot the refractive index of TiO2 and SiO2 vs lambda curve
% figure(2)
% plot(lambda,n_TiO2(:,1),'+','Color','b','Linewidth',2,'MarkerSize',20)
% hold on
% plot(lambda,n_FSi(:,1),'+','Color','r','Linewidth',2,'MarkerSize',20)
% ylabel('Magnitude')
% xline(810,'-',{'810 nm'},'LabelVerticalAlignment','middle','LabelOrientation','horizontal','Linewidth',2,'FontSize',25);
% set(gca, 'XColor','k', 'YColor','k')
% hold off
% title('Refractive index of TiO2 and SiO2')
% xlabel('Free Space Wavelength [nm]') 
% legend('n: TiO_2','n: SiO_2')
% ax = gca;
% ax.FontSize = 45;
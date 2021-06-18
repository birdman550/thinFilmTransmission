%Thin Film transmission calculator for normal incidence
%Upto 2 thin film layers on a substrate. Units of length are nanometers
%Sanket Deshpande

clear all;
close all;
digits(16)

%% determine wavelength range in nanometers
%either use a central wavelength or a range of wavelengths

%lambda=810; %central wavelength
lambda = linspace(750,850,11); %wavelength range

%% free-space refractive indices of various media
%n_Au = 0.15659+4.9908i; %gold refractive index at 810nm
n_Au = zeros(length(lambda),1);
n_Au(1) = 0.13883+4.4909i; n_Au(2) = 0.14123+4.5752i; n_Au(3) = 0.14430+4.6583i; n_Au(4) = 0.14737+4.7414i; n_Au(5) = 0.15045+4.8245i; n_Au(6) = 0.15352+4.9077i; n_Au(7) = 0.15659+4.9908i; n_Au(8) = 0.15966+5.0739i; n_Au(9) = 0.16126+5.1558i; n_Au(10) = 0.16267+5.2376i; n_Au(11) = 0.16408+5.3194i;

%n_Ti = 3.1745+4.01i; %titanium refractive index at 810 nm
n_Ti = zeros(length(lambda),1);
n_Ti(1) = 2.9838+4.0042i; n_Ti(2) = 3.0129+4.0100i; n_Ti(3) = 3.0452+4.0100i; n_Ti(4) = 3.0775+4.0100i; n_Ti(5) = 3.1098+4.0100i; n_Ti(6) = 3.1422+4.01i; n_Ti(7) = 3.1745+4.01i; n_Ti(8) = 3.2068+4.01i; n_Ti(9) = 3.2201+4.0037i; n_Ti(10) = 3.2314+3.9966i; n_Ti(11) = 3.2427+3.9896i;

%n_FSi = 1.4531; %fused silica refractive index at 810 nm
n_FSi = sqrt(1 + 0.6961663.*(lambda.^2)./(lambda.^2 - 0.0684043^2) + 0.4079426.*(lambda.^2)./(lambda.^2 - 0.1162414^2) + 0.8974794.*(lambda.^2)./(lambda.^2 - 9.896161^2) ); %Sellmeier relation from https://refractiveindex.info/?shelf=glass&book=fused_silica&page=Malitson

n_air = 1; %air refractive index

n_ITO = 1.5876+0.0059011i; %Indium Tin Oxide at 810 nm
n_ZnO = 1.9580; %Zinc oxide at 810 nm

%n_TiO2 = 2.5173; %Titanium dioxide at 810 nm
n_TiO2 = vpa(sqrt(5.913 + 0.2441./(lambda.^2 - 0.0803))); %Sellmeier relation from https://refractiveindex.info/?shelf=main&book=TiO2&page=Devore-o

%% variable declarations
n_mat = vpa(zeros(4,length(lambda))); %matrix holding refractive indices of all layers
l_mat = zeros(4,1); %matrix holding lengths of all layers. Length units are nanometers
t_mat = zeros(3,length(lambda)); %transmission coeff matrix for all interfaces
r_mat = zeros(3,length(lambda)); %reflection coeff matrix for all interfaces
phi_mat = zeros(4,length(lambda)); %phase shift matrix for all layers
T_mat = zeros(2,2,3,length(lambda)); %Transfer Matrix for interfaces
P_mat = zeros(2,2,4,length(lambda)); %Phase shift Matrix for layers
M = zeros(2,2,length(lambda)); %Total transfer matrix
for i=1:length(lambda)
    M(:,:,i) = [1 0;0 1];
end

%% defining layers in the structure
%if you intend to use less than 4 layers, you can define layers 2 and 3 to
%be made up of the same material. For example, 20 nm of Au can be defined
%as 10 nm Au in layer 1 and 10 nm Au in layer 2

%layer 1 is air with no length
n_mat(1,:) = n_air;
l_mat(1) = 0;

%layer 2 (thin film)
n_mat(2,:)=n_Au;
l_mat(2)=29;

%layer 3 (thin film)
n_mat(3,:)=n_Ti;
l_mat(3)=8;

%layer 4 (substrate)
n_mat(4,:)=n_FSi;
l_mat(4)=42;

%calculating the transfer matrix
for i=1:length(lambda)
    for a=1:3
        t_mat(a,i)= 2*n_mat(a,i)/(n_mat(a,i)+n_mat(a+1,i));
        r_mat(a,i)= (n_mat(a,i)-n_mat(a+1,i))/(n_mat(a,i)+n_mat(a+1,i));

        T_mat(:,:,a,i) = t_mat(a,i)^(-1) * [1 r_mat(a,i); r_mat(a,i) 1];
    end
end

%calculating the propagation matrix
for i=1:length(lambda)
    for b=1:3
        phi_mat(b,i) = (2*pi/lambda(i))*n_mat(b+1,i)*l_mat(b+1);
        P_mat(:,:,b,i) = [exp(-1i*phi_mat(b,i)) 0; 0 exp(1i*phi_mat(b,i))];
    end
end

%calculating the total transfer matrix
for i=1:length(lambda)
    for c=1:3
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

%plot the transmission vs lambda curve
figure(1)
plot(lambda,T.*100,'-o','Color','b','Linewidth',3)
title('Transmittance vs Free-space wavelength')
xlabel('Free Space Wavelength [nm]') 
ylabel('Transmittance [%]')
ax = gca;
ax.FontSize = 45;

%plot the phase shift vs lambda curve
figure(2)
plot(lambda,phase,'-o','Color','r','Linewidth',3)
title('Phase shift vs Free-space wavelength')
xlabel('Free Space Wavelength [nm]') 
ylabel('Phase shift [deg]')
ax = gca;
ax.FontSize = 45;


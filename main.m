%% Simplified implementation of the paper
% A. Shahmansoori, G. E. Garcia, G. Destino, G. Seco-Granados and H. Wymeersch, 
% "Position and Orientation Estimation Through Millimeter-Wave MIMO in 5G Systems," 
% in IEEE Transactions on Wireless Communications, vol. 17, no. 3, pp. 1822-1835, March 2018.
% code version 1.0, September 2018, Henk Wymeersch & Gonzalo Seco Granados, henkw@chalmers.se

close all
clear all

%% System parameters
L=4;                % number of paths (including LOS)
Rs=100;             % total BW in MHz
N=10;               % number of subcarriers 
Nt=32;              % number of TX antennas
Nr=Nt;              % number of RX antennas
Nb=Nt*2;            % number of beams in dictionary; 
Ns=20;              % number of beams sent
c=300;              % speed of light meter / us
Ts=1/Rs;            % sampling period in us
posRx=[5 1]';       % RX (user) position, TX is assumed to be in [0, 0]
alpha=0.2;            % user orientation
sigma=0.1;          % noise standard deviation


%% generate scatter points
SP=rand(L-1,2)*20-10;      % random points uniformly placed in a 20 m x 20 m area 


%% Compute Channel Parameters for L paths
TOA(1)=norm(posRx)/c;                                   % LOS TOA
AOD(1)=atan2(posRx(2),posRx(1));                        % LOS AOD
AOA(1)=atan2(posRx(2),posRx(1))-pi-alpha;               % LOS AOA
for l=1:L-1
    AOD(l+1)=atan2(SP(l,2),SP(l,1));
    AOA(l+1)=atan2(SP(l,2)-posRx(2),SP(l,1)-posRx(1))-alpha;
    TOA(l+1)=(norm(SP(l,:))+norm(posRx-SP(l,:)))/c;     % note: max distance should be below (N*Ts*c)
end

h=10*ones(1,L);                                         % some high channel gains


%% Create dictionary 
Ut=zeros(Nt,Nb);
Ur=zeros(Nr,Nb);
aa=-Nb/2:Nb/2-1;
aa=2*aa/Nb;                                 % dictionary of spatial frequencies
for m=1:Nb
    Ut(:,m)=getResponse(Nt,aa(m))*sqrt(Nt);  
    Ur(:,m)=getResponse(Nr,aa(m))*sqrt(Nr);
end


%% Generate channel: eq. (1)-(5) from the paper
H=zeros(Nr,Nt,N);
for n=1:N
    for l=1:L        
        H(:,:,n)=H(:,:,n)+h(l)*exp(-1j*2*pi*TOA(l)*(n-1)/(N*Ts))*sqrt(Nr)*getResponse(Nr,sin(AOA(l)))*sqrt(Nt)*getResponse(Nt,sin(AOD(l)))';
    end   
end


%% Visualize the beamspace channel for 1 subcarrier in AOA/AOD space
Hb=Ur'*H(:,:,1)*Ut;
mesh(asin(aa),asin(aa),abs(Hb));
xlabel('AOD'); ylabel('AOA');


%% Generate the observation and beamformers
y=zeros(Nr,Ns,N);
F=zeros(Nt,Ns,N);
for k=1:Ns    
    for n=1:N        
        F(:,k,n)=exp(1j*rand(Nt,1)*2*pi);                                           % random beamformers (note: we don't add data symbols, they are part of F)
        y(:,k,n)=H(:,:,n)*F(:,k,n)+sigma/sqrt(2)*(randn(Nr,1)+1j*randn(Nr,1));      % eq. (6) from the paper    
    end
end



%% Vectorize and generation of the basis
ybb=zeros(Nr*Ns,N);
Omega=zeros(Nr*Ns,Nb*Nb,N);
for n=1:N
    yb(:,n)=reshape(y(:,:,n),Nr*Ns,1);              % eq. (36)    
    Omega(:,:,n)=kron((Ut'*F(:,:,n)).',Ur);         % eq. (37)-(39)
end

%% run DCS-SOMP
[indices,h_hat]=DCSSOMP(yb,Omega,L);                  % the last input is the number of paths it recovers 


%% Estimate the RX position (here only with LOS)
% compute the distance
for l=1:L
    distances(l)=-mean(diff(phase(h_hat(l,:))))*(N*Ts)*c/(2*pi);                % simplified from (51)-(54)
    if (distances(l)<0)
        distances(l)=distances(l)+N*Ts*c;
    end
end

% determine LOS path index
[minval,mini]=min(distances);

% map indices back to angle
for l=1:L
    index1(l)=ceil(indices(l)/Nb);
    index2(l)=indices(l)-(index1(l)-1)*Nb;
    AOD_hat(l)=asin(aa(index1(l)));
    AOA_hat(l)=pi*sign(asin(aa(index2(l))))-asin(aa(index2(l)));
end

% compute position: 
posRx_hat = distances(mini)*[cos(AOD_hat(mini)) sin(AOD_hat(mini))]';           % eq (61)
alpha_hat= mod(AOD_hat(mini)-AOA_hat(mini)-pi,pi);
localizationError = norm(posRx_hat-posRx)                                       % in meters
orientationError = norm(alpha_hat-alpha)                                        % in rad






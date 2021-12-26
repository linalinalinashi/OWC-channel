clc;clear;close all
%----Physical/Environment Parameters----%
C=3e8; % light speed (m/s)
lx=5; ly=5; lz=3;% room dimension in metre

%----Txs Setting----%
theta=70;% semi-angle at half power
m=-log10(2)/log10(cosd(theta));% Lambertian order of emission
NumTxs=4;%Num of Txs
TPV=zeros(NumTxs,3);%initialization of Tx position metrix
TOV=zeros(NumTxs,3);%initialization of Tx orientation metrix
TPV(1,:)=[1.5 1.5 lz];% transmitter 1 position vector
TOV(1,:)=[0 0 -1];% transmitter 1 orientation vector
TPV(2,:)=[1.5 -1.5 lz];% transmitter 2 position vector
TOV(2,:)=[0 0 -1];%transmitter 2 orientation vector
TPV(3,:)=[-1.5 1.5 lz];% transmitter 3 position vector
TOV(3,:)=[0 0 -1];%transmitter 3 orientation vector
TPV(4,:)=[-1.5 -1.5 lz];% transmitter 4 position vector
TOV(4,:)=[0 0 -1];%transmitter 4 orientation vector

%----Rx Setting----%
Adet=1e-4;% detector physical area of a PD
FOV=60;% FOV of a receiver
NumRxs=1;%Num of Rxs
RPV=zeros(NumRxs,3);%initialization of receiver position metrix
ROV=zeros(NumRxs,3);%initialization of receiver orientation metrix
RPV(1,:)=[0 0 0];% receiver 1 position vector
ROV(1,:)=[0 0 1];% receiver 1 orientation vector
RPV(2,:)=[2 2 0];% receiver 2 position vector
ROV(2,:)=[0 0 1];% receiver 2 orientation vector
Nthermal=8.0185*1E-10; % thermal noise

%----Reflector Setting----%
k=1;%reflection order
rho=0.3;% reflection coefficient
Nx=lx*20; Ny=ly*20; Nz=lz*20;%reflector grid
x=-lx/2:lx/Nx:lx/2;%mesh of reflector
y=-ly/2:ly/Ny:ly/2;%mesh of reflector
z=-lz/2:lz/Nz:lz/2;%mesh of reflector
NumRefs=4;%Num of Reflector
WOV=zeros(NumRxs,3);%initialization of reflector oriantation
WPV=zeros(NumRxs,3);%initialization of reflector position
WN=zeros(NumRxs,3);%initialization of reflector mesh metrix
dA=zeros(NumRxs,1);%surface of reflector
WPV(1,:)=[-lx/2 0 0];%position of reflector 1
WOV(1,:)=[1 0 0];% wall 1 orientation vectors
WN(1,:)=[0 Ny Nz];%reflector 1 mesh metrix
dA(1,:)=lz*ly/(Nz*Ny);%reflector area
WPV(2,:)=[0 -ly/2 0];%position of wall 2
WOV(2,:)=[0 1 0];% wall 2 orientation vectors
WN(2,:)=[Nx 0 Nz];%wall 2 mesh metrix
dA(2,:)=lz*lx/(Nz*Nx);%reflector area
WPV(3,:)=[lx/2 0 0];%position of wall 3
WOV(3,:)=[-1 0 0];% wall 3 orientation vectors
WN(3,:)=[0 Ny Nz];%wall 3 mesh metrix
dA(3,:)=lz*ly/(Nz*Ny);%wall area
WPV(4,:)=[0 ly/2 0];%position of wall 4
WOV(4,:)=[0 -1 0];% wall 4 orientation vectors
WN(4,:)=[Nx 0 Nz];%wall 4 mesh metrix
dA(4,:)=lz*lx/(Nz*Nx);%reflector area

%-----Transmitted Signal-------%
S=zeros(NumTxs,16);% Initiation of signal matrix
S(1,:)=[1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0]; % emitted signal of Tx 1
S(2,:)=[0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];% emitted signal of Tx 2
S(3,:)=[0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];% emitted signal of Tx 3
S(4,:)=[0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];% emitted signal of Tx 4

%-----Dirac Function------%
Delaymax=1E-7;%dirac function duration
Delayresolution_time=1E-10;%dirac function duration
diracT=0:Delayresolution_time:Delaymax; 

%-----initialization of received signal&time----% 
Srlos=zeros(NumRxs,max(length(diracT),length(S)));
Srnlos=zeros(NumRxs,max(length(diracT),length(S)));
T=0:1:max(length(diracT),length(S))-1;
T=T*Delayresolution_time;


%-----LOS Channel impulse response---%
for rx=1:1:NumRxs
    for tx=1:1:NumTxs
        St=S(tx,:);
        D=sqrt(dot(TPV(tx,:)-RPV(rx,:),TPV(tx,:)-RPV(rx,:)));
        cosphi=lz/D;
        Delay=D/C;
        Dirac=dirac(diracT-round(Delay/Delayresolution_time)*Delayresolution_time);
        for d=1:1:length(Dirac)
            if Dirac(d)==0
                Dirac(d)=0;
            else
                Dirac(d)=1;
            end
        end
        if abs(acosd(cosphi))<=FOV
        hlos=(m+1)*Adet.*cosphi.^(m+1)./(2*pi.*D.^2)*Dirac;
        Sr_rx=filter(St,1,hlos);
        Srlos(rx,:)=Srlos(rx,:)+Sr_rx;
        end
    end
end


%-----NLOS Channel impulse response---%


for rx=1:1:NumRxs
   RPV1=RPV(rx,:);
    for tx=1:1:NumTxs
       TPV1=TPV(tx,:);
       St=S(tx,:);
        for ww=1:NumRefs
            for kk=1:WN(ww,1)+1
                for ll=1:WN(ww,2)+1
                    for mm=1:WN(ww,3)+1
                        RfPV=WPV(ww,:)+[x(kk) y(ll) z(mm)].*([1 1 1]-WOV(ww,:));
                        D1=sqrt(dot(TPV1-RfPV,TPV1-RfPV));
                        cos_phi=abs(RfPV(3)-TPV1(3))/D1;
                        cos_alpha=abs(TPV1(1)-RfPV(1))/D1;
                        D2=sqrt(dot(RfPV-RPV1,RfPV-RPV1));
                        cos_beta=abs(RfPV(1)-RPV1(1))/D2;
                        cos_psi=abs(RfPV(3)-RPV1(3))/D2;
                        tau1=(D1+D2)/C;
                        Dirac=dirac(diracT-round(tau1/Delayresolution_time)*Delayresolution_time);
                            for d=1:1:length(Dirac);
                                if Dirac(d)==0;
                                    Dirac(d)=0;
                                else
                                    Dirac(d)=1;
                                end
                            end
                            if abs(acosd(cosphi))<=FOV
                            hnlos=(m+1)*Adet*rho*dA(ww,:)*cos_phi^m*cos_alpha*cos_beta*cos_psi/(2*pi^2*D1^2*D2^2)*Dirac;
                            Sr_rx=filter(St,1,hnlos);
                            Srnlos(rx,:)=Srnlos(rx,:)+Sr_rx;
                            end
                    end
                end
            end
        end
    end
          
end

%------Received signal+noise------%
Sr_total=zeros(NumRxs,max(length(diracT),length(S)));
for rx=1:1:NumRxs
Snoise=normrnd(0,Nthermal,1,length(Srlos));
Sr_total(rx,:)=Srlos(rx,:)+Srnlos(rx,:)+Snoise;
end

%------Figure-----%
figure;
plot(T,Srnlos(1,:),'r');
figure;
plot(T,Srlos(1,:),'b');


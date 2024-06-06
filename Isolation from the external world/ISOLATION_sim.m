%% SIMULATION

% noise
sigma_p = sqrt(5/dt);
sigma_f = sqrt(5/dt);

np1 = randn(Npop,T)*sigma_p;
nf1 = randn(Npop,T)*sigma_f;
np2 = randn(Npop,T)*sigma_p;
nf2 = randn(Npop,T)*sigma_f;
npt = randn(1,T)*sigma_p;
nft = randn(1,T)*sigma_f;

inib1 = zeros(1,Npop);


%% Theta generator parameters
        
% Synaptic constants:
Ct(1,1) = 54.;  %Cep
Ct(1,2) = 54.;  %Cpe
Ct(1,3) = 54.;  %Csp
Ct(1,4) = 67.5; %Cps  
Ct(1,5) = 15.;  %Cfs
Ct(1,6) = 27.;  %Cfp
Ct(1,7) = 300.; %Cpf
Ct(1,8) = 10.;  %Cff

at(1,:)=[75*1 30*1 300]*0.85; % Reciprocal of synaptic time constants (omega) 

gaintheta = 300; % '0' for constant disinhibition 


%% L1 & L2 parameters

% Gains of synapses 
G = [5.17 4.45 57.1];

% Sigmoidal relationship:
e0 = 2.5;
r = 0.56; 
s0 = 12;

% Synaptic constants:
C = zeros(Npop,8);
C(:,1) = 54.;   % Cep
C(:,2) = 54.;   % Cpe
C(:,3) = 54.;   % Csp
C(:,4) = 67.5;  % Cps   
C(:,5) = 27.;   % Cfs
C(:,6) = 108.;  % Cfp
C(:,7) = 300.;  % Cpf
C(:,8) = 10.;   % Cff

a = ones(Npop,1)*[125 30 400];  % Reciprocal of synaptic time constants (omega) 


%% Populations inizialitation

%L1
yp1=zeros(Npop,T);
xp1=zeros(Npop,T);
vp1=zeros(Npop,1); 
zp1=zeros(Npop,T);

ye1=zeros(Npop,T);
xe1=zeros(Npop,T);
ve1=zeros(Npop,1);
ze1=zeros(Npop,T);

ys1=zeros(Npop,T);
xs1=zeros(Npop,T);
vs1=zeros(Npop,1);
zs1=zeros(Npop,T);

yf1=zeros(Npop,T);
xf1=zeros(Npop,T);
zf1=zeros(Npop,T);
vf1=zeros(Npop,1);

xl1=zeros(Npop,T);
yl1=zeros(Npop,T);

Ep1=zeros(Npop,1);
If1=zeros(Npop,1);
mf1=zeros(Npop,1);
mp1=zeros(Npop,1);

%L2
yp2=zeros(Npop,T);
xp2=zeros(Npop,T);
vp2=zeros(Npop,1); 
zp2=zeros(Npop,T);

ye2=zeros(Npop,T);
xe2=zeros(Npop,T);
ve2=zeros(Npop,1);
ze2=zeros(Npop,T);

ys2=zeros(Npop,T);
xs2=zeros(Npop,T);
vs2=zeros(Npop,1);
zs2=zeros(Npop,T);

yf2=zeros(Npop,T);
xf2=zeros(Npop,T);
zf2=zeros(Npop,T);
vf2=zeros(Npop,1);

xl2=zeros(Npop,T);
yl2=zeros(Npop,T);

Ep2=zeros(Npop,1);
If2=zeros(Npop,1);
mf2=zeros(Npop,1);
mp2=zeros(Npop,1);

%Theta generator
ypt=zeros(1,T);
xpt=zeros(1,T);
vpt=zeros(1,1); 
zpt=zeros(1,T);

yet=zeros(1,T);
xet=zeros(1,T);
vet=zeros(1,1);
zet=zeros(1,T);

yst=zeros(1,T);
xst=zeros(1,T);
vst=zeros(1,1);
zst=zeros(1,T);

yft=zeros(1,T);
xft=zeros(1,T);
zft=zeros(1,T);
vft=zeros(1,1);

xlt=zeros(1,T);
ylt=zeros(1,T);

Ept=zeros(1,1);
Ift=zeros(1,1);
mft=0;
mpt=500;


%% Isolation from the external world simulation

for k=1:T

    mp1=mediaIN(:,k);
    up1=np1(:,k)+mp1;
    uf1=nf1(:,k)+mf1;
    up2=np2(:,k)+mp2;
    uf2=nf2(:,k)+mf2;
    upt=npt(1,k)+mpt;
    uft=nft(1,k)+mft;
    
      for j = 1:Npop

        if k>1 %(i.e., D_intraL1)
        inib1(j)=Af_L1L1(j,:)*zp1(:,k-1);
        up1(j)=up1(j)+Wp_L1L1(j,:)*zp1(:,k-1);
        uf1(j)=uf1(j)+Wf_L1L1(j,:)*zp1(:,k-1);
        end
        
        if k>30 %(i,e., D_L2L1 & D_L1theta)
        up1(j)=up1(j)+Wp_L1L2(j,:)*zp2(:,k-30);
        up2(j)=up2(j)+Wp_L2L1(j,:)*zp1(:,k-30);
        up1(j)=up1(j)+gaintheta*(zpt(k-30)-5);
        end

      end
    
    % L1 post-sinaptic potentials
    vp1(:)=C(:,2).*ye1(:,k)-C(:,4).*ys1(:,k)-C(:,7).*yf1(:,k);
    ve1(:)=C(:,1).*yp1(:,k);
    vs1(:)=C(:,3).*yp1(:,k);
    vf1(:)=C(:,6).*yp1(:,k)-C(:,5).*ys1(:,k)-C(:,8).*yf1(:,k)+yl1(:,k)+inib1';
    % L1 spikes
    zp1(:,k)=2*e0./(1+exp(-r*(vp1(:)-s0))); 
    ze1(:,k)=2*e0./(1+exp(-r*(ve1(:)-s0)));
    zs1(:,k)=2*e0./(1+exp(-r*(vs1(:)-s0)));
    zf1(:,k)=2*e0./(1+exp(-r*(vf1(:)-s0)));
    % L1 new outputs
    xp1(:,k+1)=xp1(:,k)+(G(1)*a(:,1).*zp1(:,k)-2*a(:,1).*xp1(:,k)-a(:,1).*a(:,1).*yp1(:,k))*dt; 
    yp1(:,k+1)=yp1(:,k)+xp1(:,k)*dt; 
    xe1(:,k+1)=xe1(:,k)+(G(1)*a(:,1).*(ze1(:,k)+up1(:)./C(:,2))-2*a(:,1).*xe1(:,k)-a(:,1).*a(:,1).*ye1(:,k))*dt;
    ye1(:,k+1)=ye1(:,k)+xe1(:,k)*dt; 
    xs1(:,k+1)=xs1(:,k)+(G(2)*a(:,2).*zs1(:,k)-2*a(:,2).*xs1(:,k)-a(:,2).*a(:,2).*ys1(:,k))*dt;
    ys1(:,k+1)=ys1(:,k)+xs1(:,k)*dt; 
    xl1(:,k+1)=xl1(:,k)+(G(1)*a(:,1).*uf1(:)-2*a(:,1).*xl1(:,k)-a(:,1).*a(:,1).*yl1(:,k))*dt; 
    yl1(:,k+1)=yl1(:,k)+xl1(:,k)*dt; 
    xf1(:,k+1)=xf1(:,k)+(G(3)*a(:,3).*zf1(:,k)-2*a(:,3).*xf1(:,k)-a(:,3).*a(:,3).*yf1(:,k))*dt;  
    yf1(:,k+1)=yf1(:,k)+xf1(:,k)*dt; 


    % L2 post-synaptic potentials
    vp2(:)=C(:,2).*ye2(:,k)-C(:,4).*ys2(:,k)-C(:,7).*yf2(:,k);
    ve2(:)=C(:,1).*yp2(:,k);
    vs2(:)=C(:,3).*yp2(:,k);
    vf2(:)=C(:,6).*yp2(:,k)-C(:,5).*ys2(:,k)-C(:,8).*yf2(:,k)+yl2(:,k);
    % L2 spikes
    zp2(:,k)=2*e0./(1+exp(-r*(vp2(:)-s0))); 
    ze2(:,k)=2*e0./(1+exp(-r*(ve2(:)-s0)));
    zs2(:,k)=2*e0./(1+exp(-r*(vs2(:)-s0)));
    zf2(:,k)=2*e0./(1+exp(-r*(vf2(:)-s0)));
    % L2 new outputs
    xp2(:,k+1)=xp2(:,k)+(G(1)*a(:,1).*zp2(:,k)-2*a(:,1).*xp2(:,k)-a(:,1).*a(:,1).*yp2(:,k))*dt; 
    yp2(:,k+1)=yp2(:,k)+xp2(:,k)*dt; 
    xe2(:,k+1)=xe2(:,k)+(G(1)*a(:,1).*(ze2(:,k)+up2(:)./C(:,2))-2*a(:,1).*xe2(:,k)-a(:,1).*a(:,1).*ye2(:,k))*dt;
    ye2(:,k+1)=ye2(:,k)+xe2(:,k)*dt; 
    xs2(:,k+1)=xs2(:,k)+(G(2)*a(:,2).*zs2(:,k)-2*a(:,2).*xs2(:,k)-a(:,2).*a(:,2).*ys2(:,k))*dt;
    ys2(:,k+1)=ys2(:,k)+xs2(:,k)*dt; 
    xl2(:,k+1)=xl2(:,k)+(G(1)*a(:,1).*uf2(:)-2*a(:,1).*xl2(:,k)-a(:,1).*a(:,1).*yl2(:,k))*dt;
    yl2(:,k+1)=yl2(:,k)+xl2(:,k)*dt; 
    xf2(:,k+1)=xf2(:,k)+(G(3)*a(:,3).*zf2(:,k)-2*a(:,3).*xf2(:,k)-a(:,3).*a(:,3).*yf2(:,k))*dt;  
    yf2(:,k+1)=yf2(:,k)+xf2(:,k)*dt; 
    
    % Theta generator post-synaptic potentials
    vpt(:)=Ct(:,2).*yet(:,k)-Ct(:,4).*yst(:,k)-Ct(:,7).*yft(:,k); %-inib'
    vet(:)=Ct(:,1).*ypt(:,k);
    vst(:)=Ct(:,3).*ypt(:,k);
    vft(:)=Ct(:,6).*ypt(:,k)-Ct(:,5).*yst(:,k)-Ct(:,8).*yft(:,k)+ylt(:,k);
    % Theta generator spikes
    zpt(:,k)=2*e0./(1+exp(-r*(vpt(:)-s0))); 
    zet(:,k)=2*e0./(1+exp(-r*(vet(:)-s0)));
    zst(:,k)=2*e0./(1+exp(-r*(vst(:)-s0)));
    zft(:,k)=2*e0./(1+exp(-r*(vft(:)-s0)));
    % Theta generators new outputs
    xpt(:,k+1)=xpt(:,k)+(G(1)*at(1).*zpt(:,k)-2*at(1).*xpt(:,k)-at(1).*at(1).*ypt(:,k))*dt;
    ypt(:,k+1)=ypt(:,k)+xpt(:,k)*dt; 
    xet(:,k+1)=xet(:,k)+(G(1)*at(1).*(zet(:,k)+upt(:)./Ct(2))-2*at(1).*xet(:,k)-at(1).*at(1).*yet(:,k))*dt;
    yet(:,k+1)=yet(:,k)+xet(:,k)*dt; 
    xst(:,k+1)=xst(:,k)+(G(2)*at(2).*zst(:,k)-2*at(2).*xst(:,k)-at(2).*at(2).*yst(:,k))*dt;
    yst(:,k+1)=yst(:,k)+xst(:,k)*dt;
    xlt(:,k+1)=xlt(:,k)+(G(1)*at(1).*uft(:)-2*at(1).*xlt(:,k)-at(1).*at(1).*ylt(:,k))*dt; 
    ylt(:,k+1)=ylt(:,k)+xlt(:,k)*dt; 
    xft(:,k+1)=xft(:,k)+(G(3)*at(3).*zft(:,k)-2*at(3).*xft(:,k)-at(3).*at(3).*yft(:,k))*dt; 
    yft(:,k+1)=yft(:,k)+xft(:,k)*dt; 

end


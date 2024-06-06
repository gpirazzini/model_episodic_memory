%% SIMULATION

%noise:
sigma_p = sqrt(5/dt);
sigma_f = sqrt(5/dt);
%WM
np0 = randn(Npop,T)*sigma_p;
nf0 = randn(Npop,T)*sigma_f;
%L1
np1 = randn(Npop,T)*sigma_p;
nf1 = randn(Npop,T)*sigma_f;
%L2
np2 = randn(Npop,T)*sigma_p;
nf2 = randn(Npop,T)*sigma_f;
%theta
npt = randn(1,T)*sigma_p;
nft = randn(1,T)*sigma_f;
inib1 = zeros(1,Npop);


%% Working Memory parameters

% Within-column synapse time constants  <-- identical for all layers!
a = ones(Npop,1)*[125 30 400]; %in order: ae, as, af; a=1/tau (1/second)

% Gains (mV) of synapses (G)            <-- identical for all layers!
G=[5.17 4.45 57.1];  %Ge = 5.17; (for h_e)
                     %Gs = 4.45; (for h_s)
                     %Gf = 57.1; (for h_f)

% Sigmoidal relationship                <-- identical for all layers!
e0 = 2.5; % half of the maximal saturation
r = 0.56; % slope 
s0 = 12;  % center

% Synaptic constant:
C0(:,1) = 54.;   % Cep
C0(:,2) = 54.;   % Cpe
C0(:,3) = 54.;   % Csp
C0(:,4) = 67.5;  % Cps   
C0(:,5) = 27.;   % Cfs
C0(:,6) = 108.;  % Cfp
C0(:,7) = 300.;  % Cpf
C0(:,8) = 10.;   % Cff


%% Theta generator parameters 

Ct(1,1) = 54.; %Cep
Ct(1,2) = 54.; %Cpe
Ct(1,3) = 54.; %Csp
Ct(1,4) = 67.5; %Cps  
Ct(1,5) = 15.; %Cfs
Ct(1,6) = 27.; %Cfp
Ct(1,7) = 300.; %Cpf
Ct(1,8) = 10.; %Cff

at(1,:)=[75*1 30*1 300]*0.85;

gaintheta = 300;


%% L1 & L2 parameters

C = zeros(Npop,8);
C(:,1) = 54.;   % Cep
C(:,2) = 54.;   % Cpe
C(:,3) = 54.;   % Csp
C(:,4) = 67.5;  % Cps   
C(:,5) = 27.;   % Cfs
C(:,6) = 108.;  % Cfp
C(:,7) = 300.;  % Cpf
C(:,8) = 10.;   % Cff


%% Populations initialization

%WM
yp0=zeros(Npop,T);
xp0=zeros(Npop,T);
vp0=zeros(Npop,T); 
zp0=zeros(Npop,T);

ye0=zeros(Npop,T);
xe0=zeros(Npop,T);
ve0=zeros(Npop,T);
ze0=zeros(Npop,T);

ys0=zeros(Npop,T);
xs0=zeros(Npop,T);
vs0=zeros(Npop,T);
zs0=zeros(Npop,T);

yf0=zeros(Npop,T);
xf0=zeros(Npop,T);
zf0=zeros(Npop,T);
vf0=zeros(Npop,T);

xl0=zeros(Npop,T);
yl0=zeros(Npop,T);

Ee0=zeros(Npop,1);
Ep0=zeros(Npop,1);
mf0=zeros(Npop,1); 

%L1
yp1=zeros(Npop,T);
xp1=zeros(Npop,T);
vp1=zeros(Npop,T); 
zp1=zeros(Npop,T);

ye1=zeros(Npop,T);
xe1=zeros(Npop,T);
ve1=zeros(Npop,T);
ze1=zeros(Npop,T);

ys1=zeros(Npop,T);
xs1=zeros(Npop,T);
vs1=zeros(Npop,T);
zs1=zeros(Npop,T);

yf1=zeros(Npop,T);
xf1=zeros(Npop,T);
zf1=zeros(Npop,T);
vf1=zeros(Npop,T);

xl1=zeros(Npop,T);
yl1=zeros(Npop,T);

Ep1=zeros(Npop,1);
If1=zeros(Npop,1);
mf1=zeros(Npop,1);
mp1=zeros(Npop,1);

%L2
yp2=zeros(Npop,T);
xp2=zeros(Npop,T);
vp2=zeros(Npop,T); 
zp2=zeros(Npop,T);

ye2=zeros(Npop,T);
xe2=zeros(Npop,T);
ve2=zeros(Npop,T);
ze2=zeros(Npop,T);

ys2=zeros(Npop,T);
xs2=zeros(Npop,T);
vs2=zeros(Npop,T);
zs2=zeros(Npop,T);

yf2=zeros(Npop,T);
xf2=zeros(Npop,T);
zf2=zeros(Npop,T);
vf2=zeros(Npop,T);

xl2=zeros(Npop,T);
yl2=zeros(Npop,T);

Ep2=zeros(Npop,1);
If2=zeros(Npop,1);
mf2=zeros(Npop,1);
mp2=zeros(Npop,1);

%Theta generator
ypt=zeros(1,T);
xpt=zeros(1,T);
vpt=zeros(1,T); 
zpt=zeros(1,T);

yet=zeros(1,T);
xet=zeros(1,T);
vet=zeros(1,T);
zet=zeros(1,T);

yst=zeros(1,T);
xst=zeros(1,T);
vst=zeros(1,T);
zst=zeros(1,T);

yft=zeros(1,T);
xft=zeros(1,T);
zft=zeros(1,T);
vft=zeros(1,T);

xlt=zeros(1,T);
ylt=zeros(1,T);

Ept=zeros(1,1);
Ift=zeros(1,1);
mft=0;
mpt=500; %excitatory input to the Theta generator pyramidal neurons


%% Retrieval simulation 

Wp_WMWM=zeros(Npop);
IN=zeros(Npop,1);

for k=2:T-1
    %inputs to pyramidal and gaba fast (of WM):
    mp0=INPUT_WM(:,k)*5000;     
    if (sum(mp0)==0 && sum(INPUT_WM(:,k-1))>0) %if I stop receiving input...
      Wp_WMWM=diag(INPUT_WM(:,k-1))*2000; %...I keep the input in memory
    elseif  sum(mp0)~=0 %if, on the other hand, the input is nonzero...
        Wp_WMWM=zeros(Npop); %...I follow the input without self-excitation.
    end

    up0=np0(:,k)+mp0;
    uf0=nf0(:,k)+mf0;
    up1=np1(:,k)+mp1;
    uf1=nf1(:,k)+mf1;
    up2=np2(:,k)+mp2;
    uf2=nf2(:,k)+mf2;
    upt=npt(1,k)+mpt;
    uft=nft(1,k)+mft;
    
      for j = 1:Npop

        if k>D_intraL1
        up0(j)=up0(j)+Wp_WMWM(j,:)*zp0(:,k-D_intraWM);
        inib1(j)=Af_L1L1(j,:)*zp1(:,k-D_intraL1);
        up1(j)=up1(j)+Wp_L1L1(j,:)*zp1(:,k-D_intraL1);
        uf1(j)=uf1(j)+Wf_L1L1(j,:)*zp1(:,k-D_intraL1);
        end
        
        if k>D_L2L1
        up1(j)=up1(j)+Wp_L1L2(j,:)*zp2(:,k-D_L2L1);
        up2(j)=up2(j)+Wp_L2L1(j,:)*zp1(:,k-D_L1L2);
        end

        if k>D_thetaL1
        up1(j)=up1(j)+gaintheta*(zpt(k-D_thetaL1)-5);
        end

        if k>D_WML1
        up1(j)=up1(j)+Wp_L1WM(j,:)*zp0(:,k-D_WML1);
        end
        
      end
    
    %WM post-sinaptic potentials
    vp0(:,k)=C0(:,2).*ye0(:,k)-C0(:,4).*ys0(:,k)-C0(:,7).*yf0(:,k);
    ve0(:,k)=C0(:,1).*yp0(:,k);
    vs0(:,k)=C0(:,3).*yp0(:,k);
    vf0(:,k)=C0(:,6).*yp0(:,k)-C0(:,5).*ys0(:,k)-C0(:,8).*yf0(:,k)+yl0(:,k);
    %WM spikes
    zp0(:,k)=2*e0./(1+exp(-r*(vp0(:,k)-s0))); 
    ze0(:,k)=2*e0./(1+exp(-r*(ve0(:,k)-s0)));
    zs0(:,k)=2*e0./(1+exp(-r*(vs0(:,k)-s0)));
    zf0(:,k)=2*e0./(1+exp(-r*(vf0(:,k)-s0)));
    %WM new outputs
    xp0(:,k+1)=xp0(:,k)+(G(1)*a(:,1).*zp0(:,k)-2*a(:,1).*xp0(:,k)-a(:,1).*a(1).*yp0(:,k))*dt;
    yp0(:,k+1)=yp0(:,k)+xp0(:,k)*dt;
    xe0(:,k+1)=xe0(:,k)+(G(1)*a(:,1).*(ze0(:,k)+up0(:)./C0(:,2))-2*a(:,1).*xe0(:,k)-a(:,1).*a(:,1).*ye0(:,k))*dt;
    ye0(:,k+1)=ye0(:,k)+xe0(:,k)*dt;
    xs0(:,k+1)=xs0(:,k)+(G(2)*a(:,2).*zs0(:,k)-2*a(:,2).*xs0(:,k)-a(2).*a(:,2).*ys0(:,k))*dt;
    ys0(:,k+1)=ys0(:,k)+xs0(:,k)*dt;
    xl0(:,k+1)=xl0(:,k)+(G(1)*a(:,1).*uf0(:)-2*a(:,1).*xl0(:,k)-a(:,1).*a(:,1).*yl0(:,k))*dt;
    yl0(:,k+1)=yl0(:,k)+xl0(:,k)*dt;
    xf0(:,k+1)=xf0(:,k)+(G(3)*a(:,3).*zf0(:,k)-2*a(:,3).*xf0(:,k)-a(:,3).*a(:,3).*yf0(:,k))*dt;
    yf0(:,k+1)=yf0(:,k)+xf0(:,k)*dt; 

    %L1 post-sinaptic potentials
    vp1(:,k)=C(:,2).*ye1(:,k)-C(:,4).*ys1(:,k)-C(:,7).*yf1(:,k);
    ve1(:,k)=C(:,1).*yp1(:,k);
    vs1(:,k)=C(:,3).*yp1(:,k);
    vf1(:,k)=C(:,6).*yp1(:,k)-C(:,5).*ys1(:,k)-C(:,8).*yf1(:,k)+yl1(:,k)+inib1';
    %L1 spikes
    zp1(:,k)=2*e0./(1+exp(-r*(vp1(:,k)-s0))); 
    ze1(:,k)=2*e0./(1+exp(-r*(ve1(:,k)-s0)));
    zs1(:,k)=2*e0./(1+exp(-r*(vs1(:,k)-s0)));
    zf1(:,k)=2*e0./(1+exp(-r*(vf1(:,k)-s0)));
    %L1 new outputs
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

    %L2
    vp2(:,k)=C(:,2).*ye2(:,k)-C(:,4).*ys2(:,k)-C(:,7).*yf2(:,k);
    ve2(:,k)=C(:,1).*yp2(:,k);
    vs2(:,k)=C(:,3).*yp2(:,k);
    vf2(:,k)=C(:,6).*yp2(:,k)-C(:,5).*ys2(:,k)-C(:,8).*yf2(:,k)+yl2(:,k);

    zp2(:,k)=2*e0./(1+exp(-r*(vp2(:,k)-s0))); 
    ze2(:,k)=2*e0./(1+exp(-r*(ve2(:,k)-s0)));
    zs2(:,k)=2*e0./(1+exp(-r*(vs2(:,k)-s0)));
    zf2(:,k)=2*e0./(1+exp(-r*(vf2(:,k)-s0)));

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
    
    % Theta generator
    vpt(:,k)=Ct(:,2).*yet(:,k)-Ct(:,4).*yst(:,k)-Ct(:,7).*yft(:,k);
    vet(:,k)=Ct(:,1).*ypt(:,k);
    vst(:,k)=Ct(:,3).*ypt(:,k);
    vft(:,k)=Ct(:,6).*ypt(:,k)-Ct(:,5).*yst(:,k)-Ct(:,8).*yft(:,k)+ylt(:,k);

    zpt(:,k)=2*e0./(1+exp(-r*(vpt(:,k)-s0))); 
    zet(:,k)=2*e0./(1+exp(-r*(vet(:,k)-s0)));
    zst(:,k)=2*e0./(1+exp(-r*(vst(:,k)-s0)));
    zft(:,k)=2*e0./(1+exp(-r*(vft(:,k)-s0)));

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


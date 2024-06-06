%% MAIN TRAINING

tic

% parameters
dt = 0.0001;            % integration step
tend = 0.25;            % final value of the time vector
t = (0:dt:tend);         
T = length(t);          % time vector length
G = [5.17 4.45 57.1];   % synaptic gains (excitatory and pyramidal (G(1)), inhibitory slow (G(2)), and inhibitory fast (G(3))) 

% noise standard deviation
sigma_p = sqrt(5/dt);
sigma_f = sqrt(5/dt);


%% Episode list

ep = zeros(6,15);

% 1st temporal sequence)

ep1 = [1 2 3 4];
ep(1:4,1)=ep1;

ep2 = [71 70 7 8 9 10];
ep(1:6,2)=ep2;

ep3 = [11 12 13 14 15];
ep(1:5,3)=ep3;

ep4 = [16 17 18 19 20 21];
ep(1:6,4)=ep4;

ep5 = [22 23 24 25];
ep(1:4,5)=ep5;

% 2nd temporal sequence)

ep6 = [26 27 28 29];
ep(1:4,6)=ep6;

ep7 = [30 13 32 33 34 35]; 
ep(1:6,7)=ep7;

ep8 = [36 37 38 39 40];
ep(1:5,8)=ep8;

ep9 = [41 42 43 44 45 46]; 
ep(1:6,9)=ep9;

ep10 = [47 48 49 50]; 
ep(1:4,10)=ep10;

% 3rd temporal sequence)

ep11 = [51 52 53 54];
ep(1:4,11)=ep11;

ep12 = [55 13 57 58 59 60];
ep(1:6,12)=ep12;

ep13 = [61 62 44 64 65];
ep(1:5,13)=ep13;

ep14 = [66 67 68 69 70 71];
ep(1:6,14)=ep14;

ep15 = [72 73 74 75];
ep(1:4,15)=ep15;

nEp=15; %five episodes for three sequences


%% Connectivity constants between populations

C = zeros(Npop,8);
C(:,1) = 54.;   % Cep (from pyramidal to excitatory interneurons)
C(:,2) = 54.;   % Cpe
C(:,3) = 54.;   % Csp
C(:,4) = 67.5;  % Cps   
C(:,5) = 27.;   % Cfs
C(:,6) = 108.;  % Cfp
C(:,7) = 300.;  % Cpf
C(:,8) = 10.;   % Cff

a = [125 30 400]; % reciprocal of time constants (excitatory and pyramidal (a(1)), slow inhibitory (a(2)), fast inhibitory (a(3)))


%% Sigmodial relationship

e0 = 2.5; % half of the maximal saturation
r = 0.56; % slope 
s0 = 12;  % center


%% Synapses initialization

% matrices 
Wp_L1L1=zeros(Npop,Npop); % glutamatergic excitatory within L1
Wf_L1L1=zeros(Npop,Npop);  % glutamatergic inhibitory within L1
Af_L1L1=zeros(Npop,Npop);  % glutamatergic inhibitory within L1 (much faster dynamics/desynchronizing)
Wp_L1L2=zeros(Npop,Npop); % glutamatergic excitatory from L2 to L1

% max values
Wp_L1L1_max=132;
Wf_L1L1_max=15; 
Af_L1L1_max=1.5;
Wp_L1L2_max=90; 

% learning factors
gammaWp_L1L1=0.5;
gammaWf=0.02;
gammaAf=0.04;
gammaWp_L1L2=1;

% thresholds
thresh_lowWp_L1L1=0.12;
thresh_low=0.8;
thresh_up=0.6;
thresh_lowWp_L1L2=0.95;


%% Populations initialization

% L1
% pyramidal
yp1=zeros(Npop,T);
xp1=zeros(Npop,T);
vp1=zeros(Npop,1); 
zp1=zeros(Npop,T);

% excitatory
ye1=zeros(Npop,T);
xe1=zeros(Npop,T);
ve1=zeros(Npop,1);
ze1=zeros(Npop,T);

% GABAslow
ys1=zeros(Npop,T);
xs1=zeros(Npop,T);
vs1=zeros(Npop,1);
zs1=zeros(Npop,T);

% GABAfast
yf1=zeros(Npop,T);
xf1=zeros(Npop,T);
zf1=zeros(Npop,T);
vf1=zeros(Npop,1);

xl1=zeros(Npop,T);
yl1=zeros(Npop,T);


% L2
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


%% Training

for j=1:nEp %for each episode...
    
    Ep1=zeros(Npop,1);
    If1=zeros(Npop,1);
    np1 = randn(Npop,T)*sigma_p;
    nf1 = randn(Npop,T)*sigma_f;
    mp1=zeros(Npop,1);
    mf1=zeros(Npop,1);
    
    Ep2=zeros(Npop,1);
    If2=zeros(Npop,1);
    np2 = randn(Npop,T)*sigma_p;
    nf2 = randn(Npop,T)*sigma_f;
    mp2=zeros(Npop,1);
    mf2=zeros(Npop,1);

    epJ=ep(:,j);
    epJ(epJ==0)=[];
    if j~=1 && j~=6 && j~=11
        epj=ep(:,j-1);
        epj(epj==0)=[];
    end

   for k=1:T-1 %cycle over time...
        
        if j==1  || j==6 || j==11 %if it is the first episode of the sequence...

           %pyramidal and GABAfast inputs
           mp1(epJ)=4000;
           mf1(epJ)=200; 
           up1=np1(:,k)+mp1;
           uf1=nf1(:,k)+mf1;
           %(L2 receives nothing)
           up2=np2(:,k)+mp2;
           uf2=nf2(:,k)+mf2;

        else %if it is not the first episode...

           %pyramidal and GABAfast inputs
           mp1(epJ)=4000;
           mf1(epJ)=200;
           up1=np1(:,k)+mp1;
           uf1=nf1(:,k)+mf1;

           mp2(epj)=4000;
           mf2(epj)=200;
           up2=np2(:,k)+mp2;
           uf2=nf2(:,k)+mf2; 

        end

        if k>D_intraLayer
      
            Ep1=Wp_L1L1*yp1(:,k-D_intraLayer);
            If1=Wf_L1L1*yp1(:,k-D_intraLayer)+Af_L1L1*zp1(:,k-D_intraLayer); 

        end

        if k>D_extraLayer

           Ep1=Ep1+Wp_L1L2*yp2(:,k-D_extraLayer); 

        end

        %L2 post-sinaptic potentials
        vp1(:)=C(:,2).*ye1(:,k)-C(:,4).*ys1(:,k)-C(:,7).*yf1(:,k)+Ep1;
        ve1(:)=C(:,1).*yp1(:,k);
        vs1(:)=C(:,3).*yp1(:,k);
        vf1(:)=C(:,6).*yp1(:,k)-C(:,5).*ys1(:,k)-C(:,8).*yf1(:,k)+yl1(:,k)+If1;
        %L1 spikes:
        zp1(:,k)=2*e0./(1+exp(-r*(vp1(:)-s0)));
        ze1(:,k)=2*e0./(1+exp(-r*(ve1(:)-s0)));
        zs1(:,k)=2*e0./(1+exp(-r*(vs1(:)-s0)));
        zf1(:,k)=2*e0./(1+exp(-r*(vf1(:)-s0)));

        %L2 post-sinaptic potentials
        vp2(:)=C(:,2).*ye2(:,k)-C(:,4).*ys2(:,k)-C(:,7).*yf2(:,k)+Ep2;
        ve2(:)=C(:,1).*yp2(:,k);
        vs2(:)=C(:,3).*yp2(:,k);
        vf2(:)=C(:,6).*yp2(:,k)-C(:,5).*ys2(:,k)-C(:,8).*yf2(:,k)+yl2(:,k)+If2;
        %L2 spikes
        zp2(:,k)=2*e0./(1+exp(-r*(vp2(:)-s0)));
        ze2(:,k)=2*e0./(1+exp(-r*(ve2(:)-s0)));
        zs2(:,k)=2*e0./(1+exp(-r*(vs2(:)-s0)));
        zf2(:,k)=2*e0./(1+exp(-r*(vf2(:)-s0)));

        if k>300 %transient        

            %excitatory synapses
            ATT_PREw=(zp1(:,k)/(2*e0)-thresh_lowWp_L1L1)';
            ATT_PREw(ATT_PREw<0)=0;
            ATT_POSTw=(zp1(:,k)/(2*e0)-thresh_lowWp_L1L1);
            ATT_POSTw(ATT_POSTw<0)=0;
            WEIGHT_w=(Wp_L1L1_max - Wp_L1L1).*(ones(Npop,Npop)-eye(Npop));
            Wp_L1L1 = Wp_L1L1 + gammaWp_L1L1.*(ATT_POSTw * ATT_PREw) .* WEIGHT_w;

            %inhibitory synapses
            ATT_PREk=(zp1(:,k)/(2*e0) - thresh_low)';
            ATT_PREk(ATT_PREk<0)=0;
            ATT_POSTk=(zf1(:,k)/(2*e0)-thresh_low);
            ATT_POSTk(ATT_POSTk<0)=0;
            WEIGHT_k=(Wf_L1L1_max - Wf_L1L1).*(ones(Npop, Npop)-eye(Npop));
            Wf_L1L1 = Wf_L1L1 + gammaWf.*(ATT_POSTk * ATT_PREk) .* WEIGHT_k;

            %desynchronizing synapses
            ATT_PREk=(zp1(:,k)/(2*e0) - thresh_low)';
            ATT_PREk(ATT_PREk<0)=0;
            ATT_POSTa=(thresh_up-zf1(:,k)./(2*e0));
            ATT_POSTa(ATT_POSTa<0)=0;
            WEIGHT_a=(Af_L1L1_max - Af_L1L1).*(ones(Npop, Npop)-eye(Npop));
            Af_L1L1 = Af_L1L1 + gammaAf.*(ATT_POSTa * ATT_PREk) .* WEIGHT_a;

            %feedback synapses
            ATT_PRE=(zp2(:,k)/(2*e0) - thresh_lowWp_L1L2)';
            ATT_PRE(ATT_PRE<0)=0;
            ATT_POST=(zp1(:,k)/(2*e0)-thresh_lowWp_L1L2);
            ATT_POST(ATT_POST<0)=0;
            WEIGHT=(Wp_L1L2_max - Wp_L1L2).*(ones(Npop,Npop)-eye(Npop));
            Wp_L1L2 = Wp_L1L2 + gammaWp_L1L2.*(ATT_POST*ATT_PRE).*WEIGHT;

        end
        
        %new populations output:
        %L1
        %pyr
        xp1(:,k+1)=xp1(:,k)+(G(1)*a(1)*zp1(:,k)-2*a(1)*xp1(:,k)-a(1)*a(1)*yp1(:,k))*dt;
        yp1(:,k+1)=yp1(:,k)+xp1(:,k)*dt;
        %exc
        xe1(:,k+1)=xe1(:,k)+(G(1)*a(1)*(ze1(:,k)+up1(:)./C(:,2))-2*a(1)*xe1(:,k)-a(1)*a(1)*ye1(:,k))*dt;
        ye1(:,k+1)=ye1(:,k)+xe1(:,k)*dt;
        %slow
        xs1(:,k+1)=xs1(:,k)+(G(2)*a(2)*zs1(:,k)-2*a(2)*xs1(:,k)-a(2)*a(2)*ys1(:,k))*dt;
        ys1(:,k+1)=ys1(:,k)+xs1(:,k)*dt;
        %fast
        xl1(:,k+1)=xl1(:,k)+(G(1)*a(1)*uf1(:)-2*a(1)*xl1(:,k)-a(1)*a(1)*yl1(:,k))*dt;
        yl1(:,k+1)=yl1(:,k)+xl1(:,k)*dt;
        xf1(:,k+1)=xf1(:,k)+(G(3)*a(3)*zf1(:,k)-2*a(3)*xf1(:,k)-a(3)*a(3)*yf1(:,k))*dt;
        yf1(:,k+1)=yf1(:,k)+xf1(:,k)*dt;
        
        %L2
        xp2(:,k+1)=xp2(:,k)+(G(1)*a(1)*zp2(:,k)-2*a(1)*xp2(:,k)-a(1)*a(1)*yp2(:,k))*dt;
        yp2(:,k+1)=yp2(:,k)+xp2(:,k)*dt;
        xe2(:,k+1)=xe2(:,k)+(G(1)*a(1)*(ze2(:,k)+up2(:)./C(:,2))-2*a(1)*xe2(:,k)-a(1)*a(1)*ye2(:,k))*dt;
        ye2(:,k+1)=ye2(:,k)+xe2(:,k)*dt;
        xs2(:,k+1)=xs2(:,k)+(G(2)*a(2)*zs2(:,k)-2*a(2)*xs2(:,k)-a(2)*a(2)*ys2(:,k))*dt;
        ys2(:,k+1)=ys2(:,k)+xs2(:,k)*dt;
        xl2(:,k+1)=xl2(:,k)+(G(1)*a(1)*uf2(:)-2*a(1)*xl2(:,k)-a(1)*a(1)*yl2(:,k))*dt;
        yl2(:,k+1)=yl2(:,k)+xl2(:,k)*dt;
        xf2(:,k+1)=xf2(:,k)+(G(3)*a(3)*zf2(:,k)-2*a(3)*xf2(:,k)-a(3)*a(3)*yf2(:,k))*dt;
        yf2(:,k+1)=yf2(:,k)+xf2(:,k)*dt;

   end

end


%% Normalization mechanism

TOTWp_L1L1=Wp_L1L1_max*3;  %limitation of the overall neurotransmitter available
TOTWf_L1L1=Wf_L1L1_max*3;
TOTAf_L1L1=Af_L1L1_max*70;
TOTWp_L1L2=Wp_L1L2_max*4;

for i=1:Npop
    S_Wp_L1L1=sum(Wp_L1L1(i,:),2);
    if S_Wp_L1L1>TOTWp_L1L1
        Wp_L1L1(i,:)=Wp_L1L1(i,:).*(TOTWp_L1L1/S_Wp_L1L1);
    end
    S_Wf_L1L1=sum(Wf_L1L1(i,:),2); %sum of synaptic weights received from the i-th column
    if S_Wf_L1L1>TOTWf_L1L1 %only if the sum overcome the TOT...
        Wf_L1L1(i,:)=Wf_L1L1(i,:).*(TOTWf_L1L1/S_Wf_L1L1);
    end
    S_Af_L1L1=sum(Af_L1L1(i,:),2);
    if S_Af_L1L1>TOTAf_L1L1
        Af_L1L1(i,:)=Af_L1L1(i,:).*(TOTAf_L1L1/S_Af_L1L1);
    end
    S_Wp_L1L2=sum(Wp_L1L2(i,:),2);
    if S_Wp_L1L2>TOTWp_L1L2
        Wp_L1L2(i,:)=Wp_L1L2(i,:).*(TOTWp_L1L2/S_Wp_L1L2);
    end
end

toc


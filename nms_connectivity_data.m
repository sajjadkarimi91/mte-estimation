function [eeg, G, fs, tt, delays_nms] = nms_connectivity_data(T, params)

Npop = params.Npop;
Wp = params.Wp;
Wf = params.Wf;
fs = params.fs;
max_delay = params.max_delay;

G = double((Wp + Wf)>0);
G = G-diag(diag(G));
% Npop Number of ROIs

% time definition
dt=0.0001;
f_eulero=1/dt;
% tend=57; % 56 sec, because one second will be excluded for transitory effects
t=(0:dt:(100*T*dt+2));
N=length(t);



%% parameters definition
% Connectivity constants
C(:,1) = 40.*ones(1,Npop); %Cep
C(:,2) = 40.*ones(1,Npop); %Cpe
C(:,3) = 40.*ones(1,Npop); %Csp
C(:,4) = 50.*ones(1,Npop); %Cps
C(:,5) = 20.*ones(1,Npop); %Cfs
C(:,6) = 40.*ones(1,Npop); %Cfp
C(:,7) = 60.*ones(1,Npop); %Cpf
C(:,8) = 20.*ones(1,Npop); %Cff




e0 = 2.5; % Saturation value of the sigmoid
r = 0.56; % Slope of the sigmoid


step_red = round(f_eulero/fs);  % step reduction from 10000 to 256 Hz

% D=0.0166*ones(1,Npop); % Delay between regions
delays_nms = randi([1,max_delay],1,Npop);
D = dt*delays_nms.*ones(1,Npop)*step_red;% Delay between regions

a=[75 30 300 ]; % Reciprocal of synaptic time constants (\omega)

SG=[5.17 4.45 57.1]; % Synaptic gains

%% Simulation


sigma = sqrt(9/dt); % Standard deviation of the input noise
np = randn(Npop,N)*sigma; % input noise to excitatory neurons
nf = randn(Npop,N)*sigma; % input noise to inhibitory neurons

% defining equations of a single ROI
yp=zeros(Npop,N);
xp=zeros(Npop,N);
vp=zeros(Npop,1);
zp=zeros(Npop,N);
ye=zeros(Npop,N);
xe=zeros(Npop,N);
ve=zeros(Npop,1);
ze=zeros(Npop,N);
ys=zeros(Npop,N);
xs=zeros(Npop,N);
vs=zeros(Npop,1);
zs=zeros(Npop,N);
yf=zeros(Npop,N);
xf=zeros(Npop,N);
zf=zeros(Npop,N);
vf=zeros(Npop,1);
xl=zeros(Npop,N);
yl=zeros(Npop,N);



m = zeros(Npop,1); % mean value of the input noise

kmax=round(max(D)/dt);

for k=1:N-1
    up=np(:,k)+m; % input of exogenous contributions to excitatory neurons
    uf=nf(:,k);  % input of exogenous contributions to inhibitory neurons
    
    if(k>kmax)
        for i=1:Npop
            up(i)=up(i)+Wp(:,i)'*zp(:,round(k-D(i)/dt));
            uf(i)=uf(i)+Wf(:,i)'*zp(:,round(k-D(i)/dt));
        end
    end
    
    % post-synaptic membrane potentials
    vp(:)=C(:,2).*ye(:,k)-C(:,4).*ys(:,k)-C(:,7).*yf(:,k);
    ve(:)=C(:,1).*yp(:,k);
    vs(:)=C(:,3).*yp(:,k);
    vf(:)=C(:,6).*yp(:,k)-C(:,5).*ys(:,k)-C(:,8).*yf(:,k)+yl(:,k);
    
    % average spike density
    zp(:,k)=2*e0./(1+exp(-r*(vp(:))))-e0;
    ze(:,k)=2*e0./(1+exp(-r*(ve(:))))-e0;
    zs(:,k)=2*e0./(1+exp(-r*(vs(:))))-e0;
    zf(:,k)=2*e0./(1+exp(-r*(vf(:))))-e0;
    
    % post synaptic potential for pyramidal neurons
    xp(:,k+1)=xp(:,k)+(SG(1)*a(1)*zp(:,k)-2*a(1)*xp(:,k)-a(1)*a(1)*yp(:,k))*dt;
    yp(:,k+1)=yp(:,k)+xp(:,k)*dt;
    
    % post synaptic potential for excitatory interneurons
    xe(:,k+1)=xe(:,k)+(SG(1)*a(1)*(ze(:,k)+up(:)./C(:,2))-2*a(1)*xe(:,k)-a(1)*a(1)*ye(:,k))*dt;
    ye(:,k+1)=ye(:,k)+xe(:,k)*dt;
    
    % post synaptic potential for slow inhibitory interneurons
    xs(:,k+1)=xs(:,k)+(SG(2)*a(2)*zs(:,k)-2*a(2)*xs(:,k)-a(2)*a(2)*ys(:,k))*dt;
    ys(:,k+1)=ys(:,k)+xs(:,k)*dt;
    
    % post synaptic potential for fast inhibitory interneurons
    xl(:,k+1)=xl(:,k)+(SG(1)*a(1)*uf(:)-2*a(1)*xl(:,k)-a(1)*a(1)*yl(:,k))*dt;
    yl(:,k+1)=yl(:,k)+xl(:,k)*dt;
    xf(:,k+1)=xf(:,k)+(SG(3)*a(3)*zf(:,k)-2*a(3)*xf(:,k)-a(3)*a(3)*yf(:,k))*dt;
    yf(:,k+1)=yf(:,k)+xf(:,k)*dt;
    
    
end

% 3 ROIs data generation
start = 10000; % exclusion of the first second due to a possible transitory
eeg = diag(C(:,2))*ye(:,start:step_red:end)-diag(C(:,4))*ys(:,start:step_red:end)-diag(C(:,7))*yf(:,start:step_red:end);
eeg = eeg(:,1:T);

% eeg = diag(C(:,2))*ye(:,start:end)-diag(C(:,4))*ys(:,start:end)-diag(C(:,7))*yf(:,start:end);

tt=t(start:step_red:end); % time vector


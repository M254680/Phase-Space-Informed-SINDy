clear all, clc, close all
set(0, 'DefaultAxesFontSize', 14);   % Set default font size for plots
set(0, 'DefaultTextFontSize', 14);   % Set default font size for text

% Select system: 
% eq = 1 → Two-Attractor, 2 → Duffing oscillator, 3 → Lorenz system
eq = 1;  



% Parameters for the  system
    % Lorenz's parameters (chaotic)
        sigma = 10;  
        beta = 8/3;
        rho = 28;

    % Duffing system parameters
        lambda_D=1.3;
        gamma=.5;


% Time settings
dt = 0.001;  % Time step
t_final = 25;  % Total time
t_span = 0.001:dt:t_final;  % Time vector

% Initial conditions
    %Two Attractor
        x01=[-2;-2];  
        % x02=[-1/2;-5];
        %   x03=[1/2;2];
          % x03=[3;-4];

          %Test case the breaks 
% x01=[3;-4];


       x05=[5;1/25];
       % x05=[2;2];
        
     %Duffing

     
%IC 1 and 2 on different lobs
     % x01=[-1;7/4];  
        % x02=[-3/2;-5];
          % x03=[0;-11/3];
          % x04=[3;-4];
%2IC with 1 lobe
     % x01=[3;-4];
     %  x02=[-3/2;-5];

%Tight on left
% x01=[-2.3;-0];

%Tight on Right
 % x01=[0.81;.75];

       % x05=[-5;2/5];  
     %Lorenz
        % x01=[-2;-2;0];  %right lobe

        %  x01=[2;-5;7];
        % x02=[1/2;2;-2];  %left lobe
        % x03=[-1/2;-5;3];
        % x04=[3;-4;1];
        % 
        % x05=[2;5;15];

x01=[3,3];
x02=[-2,-1/2];
x03=[4.8,-1.001];

x05=[-4,1/2];

% Integrate the True System
n=2; %number of state variables
options = odeset('RelTol',1e-12,'AbsTol',1e-12*ones(1,n));
[tT1,xT1]=ode45(@(t,x)System(t,x, lambda_D, gamma, sigma, beta, rho,eq),t_span,x01,options);

%Plot the Phase Portrait
figure
hold on
title('Trajectories of Two-Attractor')
plot(xT1(:,1),xT1(:,2))
grid on
legend('x01','x02','x03','IC4')
xlim([-5/2,5/2])
ylim([-5/2,5/2])

%% Add Noise to Trajectory

rng(7401) % Seed for reproducibility
%Noise method "+- <= percentage of noise.
        noise_level = 0.005;  % 0.5% noise

        % % Generate random noise in the range [-noise_level*100%, +noise_level*100 %] of each element in xT1
        % %take value between 0 &1, subtracts .5 and multiples by 2 to get +-
        noise1 = (rand(size(xT1)) - 0.5) * 2 * noise_level .* xT1;
        

        % Add the noise to the original matrix
        xT1N = xT1 + noise1;   % Contaminated trajectory
 
    % Compare clean vs noisy
           figure
          hold on
        plot(tT1,xT1N(:,1))
        plot(tT1,xT1(:,1))
        title("Contaminated Data vs Uncontaminated Data")
        xlabel('t')
        ylabel('x')
        xlim([1,4])
        % ylim([-1.25,-.75])
        legend('Contaminated','Uncontaminated')
%% Set SINDy Library Parameters
%libary parameters
polyorder = 3; % Polynomial order of library (up to cubic terms)
usesine = 0; % Do not use sine terms

%dimension set
n=n; %preset already for ODE45


 %Defing True Coefficents
%  XiT=SystemXiT( lambda_D, gamma, sigma, beta, rho,eq);
% poolDataLIST({'x','y'},XiT,n,polyorder,usesine);
% poolDataLIST({'x','y','z'},XiT,n,polyorder,usesine);


%% Savitzky Golay Derivative
% Set Savitzky-Golay filter parameters
windowSize = 25;   % Must be odd
polynomialOrder = 3;  % Polynomial order, e.g., cubic; dicates how big of a derivate you can take. this can take up to 3rd derivative

% Smooth the signal
% smoothedSignal = sgolayfilt(xT1N, polynomialOrder, windowSize);
% plot(tT,smoothedSignal(:,1),'LineWidth',1.5)
% xlim([6,8])
%% Determining Savitzky-Golay Derivative Estimation Coefficents
[b, g] = sgolay(polynomialOrder, windowSize);

%Estimate derivatives using convolution
dxT1N=zeros(size(xT1N));
for i=1:n
    dxT1N(:,i)=conv(xT1N(:,i), factorial(1)*g(:,2)/(-dt)^1, 'same');
end

% Truncate boundary points due to filtering edge effects
dxT1N_SG=dxT1N((windowSize+1)/2:end-((windowSize-1)/2),:);


%% Finite Difference Derivative Estimation

for i=1:length(xT1N)-1
    dxT1N_DD(i,:)=(xT1N(i+1,:)-xT1N(i,:))/(dt); %suppose time between snapshots is 0.001 second
end

%% Align Data with Derivative Sizes
 
 xT1N_SG=xT1N((windowSize+1)/2:end-((windowSize-1)/2),:);

 xT1N_DD=xT1N(1:length(xT1N)-1,:);

%% Build the Library of Candidate Terms
Theta_SG= poolData(xT1N_SG,n,polyorder,usesine);
Theta_DD= poolData(xT1N_DD,n,polyorder,usesine);


%% Sparse Regression to Identify Dynamics
lambda = 0.25;      % Sparsification threshold

iters=25;
Xi_SG = sparsifyDynamics(Theta_SG,dxT1N_SG,lambda,n,iters);

Xi_DD = sparsifyDynamics(Theta_DD,dxT1N_DD,lambda,n,iters);
% Display the identified models
poolDataLIST({'x','y'},Xi_SG,n,polyorder,usesine);
% poolDataLIST({'x','y','z'},XiC,n,polyorder,usesine);

poolDataLIST({'x','y'},Xi_DD,n,polyorder,usesine);
% poolDataLIST({'x','y','z'},Xi1,n,polyorder,usesine);

% %% Results
% 
% [tC,xi_SG_1]=ode45(@(t,x)sparseGalerkin(t,x,Xi_SG,polyorder,usesine),t_span,x01,options); 
%  [t1,xi_DD_1]=ode45(@(t,x)sparseGalerkin(t,x,Xi_DD,polyorder,usesine),t_span,x01,options); 
% 
% 
% 
% figure
% hold on
% title('Two-Attractor X-Time Series of True and SINDys on IC1 (Training Validation)')
% plot(tT1,xT1(:,1),'k','LineWidth',1.75)
% plot(tC,xi_SG_1(:,1),'r--','LineWidth',1.5)
% plot(t1,xi_DD_1(:,1),'b--','LineWidth',1.25)
% 
% legend('True','SG','DD')
% ylabel('x')
% xlabel('t')
% 
% RRMSE_XiC1_1=RRMSE(xT1, xi_SG_1);
% RRMSE_XiC2_1=RRMSE(xT1, xi_DD_1);
% RRMSEall_1=[RRMSE_XiC1_1, RRMSE_XiC2_1]
% 
% 
% [tT5,xT5]=ode45(@(t,x)System(t,x, lambda_D, gamma, sigma, beta, rho,eq),t_span,x05,options);
% [tC,xi_SG_5]=ode45(@(t,x)sparseGalerkin(t,x,Xi_SG,polyorder,usesine),t_span,x05,options); 
% [t1,xi_DD_5]=ode45(@(t,x)sparseGalerkin(t,x,Xi_DD,polyorder,usesine),t_span,x05,options); 
% 
%  figure
% hold on
% title('Two-Attractor X-Time Series of True and SINDys on IC5 (Testing Validation)')
% plot(tT5,xT5(:,1),'k','LineWidth',1.75)
% plot(tC,xi_SG_5(:,1),'r--','LineWidth',1.5)
% plot(t1,xi_DD_5(:,1),'b--','LineWidth',1.25)
% 
% legend('True','SG','DD')
% ylabel('x')
% xlabel('t')
% 
% RRMSE_XiC1_5=RRMSE(xT5, xi_SG_5)
% RRMSE_XiC2_5=RRMSE(xT5, xi_DD_5);
% RRMSEall_5=[RRMSE_XiC1_5, RRMSE_XiC2_5]

%%
function Error=RRMSE(xT, xP)
Error=rmse(xT,xP)/vecnorm(xT);
end

function SparsityValue=Sparsity(XiT,XiP)
PercentDifference = zeros(size(XiT));
PercentDifference(XiT == 0 & XiP ~= 0) = 1;
PercentDifference(XiT ~= 0 & XiP == 0) = 1;

mask = (XiT ~= 0 & XiP ~= 0);
percentDiff = abs(XiT - XiP) ./abs(XiT); % Element-wise percent difference
PercentDifference(mask) = min(percentDiff(mask), 1); % Cap at 1
SparsityValue=sum(PercentDifference,'all')/(numel(XiT));

end

function XiT=SystemXiT( lambda_D, gamma, sigma, beta, rho,eq)

switch eq 
    case 1
        XiT =[ 
            0 0;
            1 0;
            0 -1;
            0 0;
            0 0;
            0 0;
            -1 0;
            0 0;
            0 0;
            0 0
        ];
    case 2
        XiT=[
            0 0;
            0 1;
            1 -gamma;
            0 -lambda_D;
            0 0;
            0 0;
            0 -1;
            0 0;
            0 0;
            0 0
            ];
    case 3
%poly order 2
XiT=[
             0 0 0;
            -sigma rho 0;
            sigma -1 0;
            0 0 -beta;
            0 0 0;
            0 -1 1;
            0 0 0;
            0 0 0;
            0 0 0;
            0 0 0;]

        %poly Order 3
        % XiT=[
        %      0 0 0;
        %     -sigma rho 0;
        %     sigma -1 0;
        %     0 0 -beta;
        %     0 0 0;
        %     0 -1 1;
        %     0 0 0;
        %     0 0 0;
        %     0 0 0;
        %     0 0 0;
        %     0 0 0;
        %     0 0 0;
        %     0 0 0;
        %     0 0 0;
        %     0 0 0;
        %     0 0 0;
        %     0 0 0;
        %     0 0 0;
        %     0 0 0;
        %     0 0 0
        %     ];
end

end

function dy = System(t,y, lambda_D, gamma, sigma, beta, rho,eq)
switch eq

    case 1
        dy = [ y(1)-y(1).^3;
    -y(2)
];

    case 2
dy = [ y(2);
    -y(1)*(-1+lambda_D*y(1)+y(1).^2)-gamma*y(2)
];
    case 3
dy = [
sigma*(y(2)-y(1));
y(1)*(rho-y(3))-y(2);
y(1)*y(2)-beta*y(3);
];
    case 4
dy = [-y(1);
    y(2).^2*(1-y(2).^2)
    ];
end


end



function dy = sparseGalerkin(t,y,ahat,polyorder,usesine)
% Copyright 2015, All Rights Reserved
% Code by Steven L. Brunton
% For Paper, "Discovering Governing Equations from Data: 
%        Sparse Identification of Nonlinear Dynamical Systems"
% by S. L. Brunton, J. L. Proctor, and J. N. Kutz

yPool = poolData(y',length(y),polyorder,usesine);
dy = (yPool*ahat)';
end



function yout = poolData(yin,nVars,polyorder,usesine)
% Copyright 2015, All Rights Reserved
% Code by Steven L. Brunton
% For Paper, "Discovering Governing Equations from Data: 
%        Sparse Identification of Nonlinear Dynamical Systems"
% by S. L. Brunton, J. L. Proctor, and J. N. Kutz

n = size(yin,1);
% yout = zeros(n,1+nVars+(nVars*(nVars+1)/2)+(nVars*(nVars+1)*(nVars+2)/(2*3))+11);

ind = 1;
% poly order 0
yout(:,ind) = ones(n,1);
ind = ind+1;

% poly order 1
for i=1:nVars
    yout(:,ind) = yin(:,i);
    ind = ind+1;
end

if(polyorder>=2)
    % poly order 2
    for i=1:nVars
        for j=i:nVars
            yout(:,ind) = yin(:,i).*yin(:,j);
            ind = ind+1;
        end
    end
end

if(polyorder>=3)
    % poly order 3
    for i=1:nVars
        for j=i:nVars
            for k=j:nVars
                yout(:,ind) = yin(:,i).*yin(:,j).*yin(:,k);
                ind = ind+1;
            end
        end
    end
end

if(polyorder>=4)
    % poly order 4
    for i=1:nVars
        for j=i:nVars
            for k=j:nVars
                for l=k:nVars
                    yout(:,ind) = yin(:,i).*yin(:,j).*yin(:,k).*yin(:,l);
                    ind = ind+1;
                end
            end
        end
    end
end

if(polyorder>=5)
    % poly order 5
    for i=1:nVars
        for j=i:nVars
            for k=j:nVars
                for l=k:nVars
                    for m=l:nVars
                        yout(:,ind) = yin(:,i).*yin(:,j).*yin(:,k).*yin(:,l).*yin(:,m);
                        ind = ind+1;
                    end
                end
            end
        end
    end
end

if(usesine)
    for k=1:10;
        yout = [yout sin(k*yin) cos(k*yin)];
    end
end
end

function Xi = sparsifyDynamics(Theta,dXdt,lambda,n,runs)
% Copyright 2015, All Rights Reserved
% Code by Steven L. Brunton
% For Paper, "Discovering Governing Equations from Data: 
%        Sparse Identification of Nonlinear Dynamical Systems"
% by S. L. Brunton, J. L. Proctor, and J. N. Kutz

% compute Sparse regression: sequential least squares
Xi = Theta\dXdt;  % initial guess: Least-squares

% lambda is our sparsification knob.
for k=1:runs
    smallinds = (abs(Xi)<lambda);   % find small coefficients
    Xi(smallinds)=0;                % and threshold
    for ind = 1:n                   % n is state dimension
        biginds = ~smallinds(:,ind);
        % Regress dynamics onto remaining terms to find sparse Xi
        Xi(biginds,ind) = Theta(:,biginds)\dXdt(:,ind); 
    end
end
end

function yout = poolDataLIST(yin,ahat,nVars,polyorder,usesine)
% Copyright 2015, All Rights Reserved
% Code by Steven L. Brunton
% For Paper, "Discovering Governing Equations from Data: 
%        Sparse Identification of Nonlinear Dynamical Systems"
% by S. L. Brunton, J. L. Proctor, and J. N. Kutz

n = size(yin,1);

ind = 1;
% poly order 0
yout{ind,1} = ['1'];
ind = ind+1;

% poly order 1
for i=1:nVars
    yout(ind,1) = yin(i);
    ind = ind+1;
end

if(polyorder>=2)
    % poly order 2
    for i=1:nVars
        for j=i:nVars
            yout{ind,1} = [yin{i},yin{j}];
            ind = ind+1;
        end
    end
end

if(polyorder>=3)
    % poly order 3
    for i=1:nVars
        for j=i:nVars
            for k=j:nVars
                yout{ind,1} = [yin{i},yin{j},yin{k}];
                ind = ind+1;
            end
        end
    end
end

if(polyorder>=4)
    % poly order 4
    for i=1:nVars
        for j=i:nVars
            for k=j:nVars
                for l=k:nVars
                    yout{ind,1} = [yin{i},yin{j},yin{k},yin{l}];
                    ind = ind+1;
                end
            end
        end
    end
end

if(polyorder>=5)
    % poly order 5
    for i=1:nVars
        for j=i:nVars
            for k=j:nVars
                for l=k:nVars
                    for m=l:nVars
                        yout{ind,1} = [yin{i},yin{j},yin{k},yin{l},yin{m}];
                        ind = ind+1;
                    end
                end
            end
        end
    end
end

if(usesine)
    for k=1:10;
        yout{ind,1} = ['sin(',num2str(k),'*yin)'];
        ind = ind + 1;
        yout{ind,1} = ['cos(',num2str(k),'*yin)'];
        ind = ind + 1;
    end
end


output = yout;
newout(1) = {''};
for k=1:length(yin)
    newout{1,1+k} = [yin{k},'dot'];
end
% newout = {'','xdot','ydot','udot'};
for k=1:size(ahat,1)
    newout(k+1,1) = output(k);
    for j=1:length(yin)
        newout{k+1,1+j} = ahat(k,j);
    end
end
newout
end
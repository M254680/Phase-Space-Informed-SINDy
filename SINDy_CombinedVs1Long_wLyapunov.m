%% Initialization and System Selection
clear all, clc, close all
set(0, 'DefaultAxesFontSize', 14); 
set(0, 'DefaultTextFontSize', 14); 

% Select the system: 1 = Two-Attractor, 2 = Duffing, 3 = Lorenz
eq=1;

if eq==3
    n=3; %number of observables
else
    n=2; %number of observables
end

%% Define System Parameters
    % Lorenz's parameters (chaotic)
        sigma = 10;  
        beta = 8/3;
        rho = 28;

    %Duffing Paraters
        lambda_D=1.3;
        gamma=.5;


%% Time Configuration
dt = 0.001;  % Time step
t_final = 10;  % Total time
t_span = 0.001:dt:t_final;  % Time vector
t_span_L = 0.001:dt:t_final*2;  % Time vector

%% Define Initial Conditions
            %Two Attractor
            if eq==1
 x01=    [-2;-2];            %L
 x02=     [-0.5,-5] ;          %L
 x03=    [3;-4]   ;   %R

 x05=    [5;1/25]   ;   %R
             %Duffing
            elseif eq==2
 x01=    [-2.5;0];            %L
 x02=     [-1.5,-5] ;          %L
 x03=    [2.4;2]   ;   %R

 x05=    [-2;2.5]   ;   %R

             %Lorenz
             else
 x01=    [5;15;20];            %I
  x02=     [-1.5;-5;5] ;          %I
  x03=    [50;60;80]   ;   %O

  x05=    [40;30;20]   ;   %O
  end


% Integrate the system
options = odeset('RelTol',1e-12,'AbsTol',1e-12*ones(1,n));
[tT1_L,xT1_Long]=ode45(@(t,x)System(t,x, lambda_D, gamma, sigma, beta, rho,eq),t_span_L,x01,options);
[tT1,xT1]=ode45(@(t,x)System(t,x, lambda_D, gamma, sigma, beta, rho,eq),t_span,x01,options);
[tT2,xT2]=ode45(@(t,x)System(t,x, lambda_D, gamma, sigma, beta, rho,eq),t_span,x02,options);
 [tT3,xT3]=ode45(@(t,x)System(t,x, lambda_D, gamma, sigma, beta, rho,eq),t_span,x03,options);

[tT5,xT5]=ode45(@(t,x)System(t,x, lambda_D, gamma, sigma, beta, rho,eq),t_span,x05,options);

% Plot Phase Portrait
figure
hold on
if eq==3
plot3(xT1(:,1),xT1(:,2),xT1(:,3),'LineWidth',1.75)
plot3(xT2(:,1),xT2(:,2),xT2(:,3),'LineWidth',1.75)
plot3(xT3(:,1),xT3(:,2),xT3(:,3),'LineWidth',1.75)
plot3(xT5(:,1),xT5(:,2),xT5(:,3),'LineWidth',1.75)
 xlabel('x')
        ylabel('y')
        zlabel('z')
view(-131,8)
else
plot(xT1(:,1),xT1(:,2),'LineWidth',1.75)
plot(xT2(:,1),xT2(:,2),'LineWidth',1.75)
plot(xT3(:,1),xT3(:,2),'LineWidth',1.75)
plot(xT5(:,1),xT5(:,2),'LineWidth',1.75)
xline(0, 'Color', 'k', 'LineWidth', 1); % Draw line for Y axis.
yline(0, 'Color', 'k', 'LineWidth', 1); % Draw line for X axis.
    xlabel('x')
        ylabel('y')
end


if eq==1
xlim([-2,2]) %For Two Attractor
ylim([-2,2])%For Two Attractor
    xlabel('x')
        ylabel('y')
end
if eq==2;
xlim([-4,4]) %For Duffing
ylim([-4,4])%For Duffing
    xlabel('x')
        ylabel('y')
end

grid on
% title('Trajectories of IC1, IC2, IC3, and IC5')
% legend('IC1','IC2', 'IC3', 'IC5')
%% Lyapunov Exponent
% dim = 3;
% [~,lag] = phaseSpaceReconstruction(xT2(:,1),[],dim);
% eRange=200;
% lyapunovExponent(xT1(:,1),1/dt,dim, 'ExpansionRange',eRange)
% Kmin = 220;
% Kmax = 400;
% lyapExp = lyapunovExponent(xT1(:,1),1/dt,lag,dim,'ExpansionRange',[Kmin Kmax])

%% Add Noise to Trajectories

rng(7401) %seed for replication
%Noise method 1 "+- <= percentage of noise.
        noise_level = 0.05; % 5% noise

%calculate noise
        noise1 = (rand(size(xT1)) - 0.5) * 2 * noise_level .* xT1;
        noise2 = (rand(size(xT2)) - 0.5) * 2 * noise_level .* xT2;
        noise3 = (rand(size(xT3)) - 0.5) * 2 * noise_level .* xT3;
        
         noise1_L = (rand(size(xT1_Long)) - 0.5) * 2 * noise_level .* xT1_Long;

        % Add the noise to the original matrix
        xT1N = xT1 + noise1;
        xT2N=xT2+noise2;
        xT3N=xT3+noise3;

        xT1N_Long = xT1_Long + noise1_L;
    
        %plot noise and normal
           figure
          hold on
        plot(tT1,xT1N(:,1))
        plot(tT1,xT1(:,1))
        title("Contaminated Data vs Uncontaminated Data")
        xlabel('t')
        ylabel('x')
        xlim([5,6])
        ylim([-1.25,-.75])
        legend('Contaminated','Uncontaminated')
%% generate Library SINDy Parameters
%libary parameters
if eq==3
polyorder = 2;
else
polyorder = 3;
end
usesine = 0;

%dimension set
n=n; %preset already for ODE45


%% Define True Coefficients for Lorenz System
 XiT=SystemXiT( lambda_D, gamma, sigma, beta, rho,eq);
 if n==2
poolDataLIST({'x','y'},XiT,n,polyorder,usesine);
 else
poolDataLIST({'x','y','z'},XiT,n,polyorder,usesine);
 end

%% Savitzky-Golay Derivative Estimation
% Set Savitzky-Golay filter parameters
windowSize = 25;  % Odd number, e.g., 11 points in the window
polynomialOrder = 3;  % Polynomial order, e.g., cubic; dicates how big of a derivate you can take. this can take up to 3rd derivative

% Smooth the signal
% smoothedSignal = sgolayfilt(xT1N, polynomialOrder, windowSize);
% plot(tT,smoothedSignal(:,1),'LineWidth',1.5)
% xlim([6,8])
% Savitzky-Golay differentiation coefficients
[b, g] = sgolay(polynomialOrder, windowSize);

% Compute derivatives for each noisy trajectory
dxT1N=zeros(size(xT1N));
for i=1:n
    dxT1N(:,i)=conv(xT1N(:,i), factorial(1)*g(:,2)/(-dt)^1, 'same');
end
dxT1N=dxT1N((windowSize+1)/2:end-((windowSize-1)/2),:);

dxT2N=zeros(size(xT2N));
for i=1:n
    dxT2N(:,i)=conv(xT2N(:,i), factorial(1)*g(:,2)/(-dt)^1, 'same');
end
dxT2N=dxT2N((windowSize+1)/2:end-((windowSize-1)/2),:);

dxT3N=zeros(size(xT3N));
for i=1:n
    dxT3N(:,i)=conv(xT3N(:,i), factorial(1)*g(:,2)/(-dt)^1, 'same');
end
dxT3N=dxT3N((windowSize+1)/2:end-((windowSize-1)/2),:);


dxT1N_Long=zeros(size(xT1N_Long));
for i=1:n
    dxT1N_Long(:,i)=conv(xT1N_Long(:,i), factorial(1)*g(:,2)/(-dt)^1, 'same');
end
dxT1N_Long=dxT1N_Long((windowSize+1)/2:end-((windowSize-1)/2),:);

%% Combining SG 
 % Trim unfiltered fringe
 xT1N=xT1N((windowSize+1)/2:end-((windowSize-1)/2),:);
 xT2N=xT2N((windowSize+1)/2:end-((windowSize-1)/2),:);
 xT3N=xT3N((windowSize+1)/2:end-((windowSize-1)/2),:);

 xT1N_Long=xT1N_Long((windowSize+1)/2:end-((windowSize-1)/2),:);
%% Combine State and Derivative Data for Library Construction
 xTC_LL=[xT1N; xT2N];
dxTC_LL=[dxT1N; dxT2N];

 xTC_LR=[xT1N; xT3N];
dxTC_LR=[dxT1N; dxT3N];

%% Build Library of Candidate Functions
Theta_LL= poolData(xTC_LL,n,polyorder,usesine);
Theta_LR= poolData(xTC_LR,n,polyorder,usesine);
Theta_Long= poolData(xT1N_Long,n,polyorder,usesine);


%% Sparse Regression to Identify Dynamics
lambda = .25;      % sparsity coefficent
iters=25;
Xi_LL = sparsifyDynamics(Theta_LL,dxTC_LL,lambda,n,iters);
 if n==2
 poolDataLIST({'x','y'},Xi_LL,n,polyorder,usesine);
 else poolDataLIST({'x','y','z'},Xi_LL,n,polyorder,usesine); end

 Xi_LR = sparsifyDynamics(Theta_LR,dxTC_LR,lambda,n,iters);
 if n==2
 poolDataLIST({'x','y'},Xi_LR,n,polyorder,usesine);
 else poolDataLIST({'x','y','z'},Xi_LR,n,polyorder,usesine); end

Xi1_Long = sparsifyDynamics(Theta_Long,dxT1N_Long,lambda,n,iters);
 if n==2
 poolDataLIST({'x','y'},Xi1_Long,n,polyorder,usesine);
 else poolDataLIST({'x','y','z'},Xi1_Long,n,polyorder,usesine); end
%% Evaluate and Plot Results on Training IC
[tC_LL,xC_LL]=ode45(@(t,x)sparseGalerkin(t,x,Xi_LL,polyorder,usesine),t_span,x01,options); 
[tC_LR,xC_LR]=ode45(@(t,x)sparseGalerkin(t,x,Xi_LR,polyorder,usesine),t_span,x01,options); 
 [t1,x1_Long]=ode45(@(t,x)sparseGalerkin(t,x,Xi1_Long,polyorder,usesine),t_span,x01,options); 

Sparsity_XiC_LL=Sparsity(XiT,Xi_LL)
Sparsity_XiC_LR=Sparsity(XiT,Xi_LR)
Sparsity_XiC_Long=Sparsity(XiT,Xi1_Long)
% Sparsity_all=[Sparsity_XiC_LL,Sparsity_XiC_LR, Sparsity_XiC_Long]

figure
hold on
plot(tT1,xT1(:,1),'k','LineWidth',1.95)
plot(tC_LL,xC_LL(:,1),'r--','LineWidth',1.75)
plot(tC_LR,xC_LR(:,1),'g--','LineWidth',1.5)
plot(t1,x1_Long(:,1),'b--','LineWidth',1.25)
if eq==1
title('Two-Attractor X-Time Series of True and SINDys on IC1 (Training Validation)')
legend('True','LL', 'LR','One Long')
end
if eq==2
title('Duffing X-Time Series of True and SINDys on IC1 (Training Validation)')
legend('True','LL', 'LR','One Long')
end
if eq==3
title('Lorenz X-Time Series of True and SINDys on IC1 (Training Validation)')
legend('True','II', 'IO','I')
end



ylabel('x')
xlabel('t')

RRMSE_Xi_LL_1=RRMSE(xT1, xC_LL)
RRMSE_Xi_LR_1=RRMSE(xT1, xC_LR)
RRMSE_Xi_Long_1=RRMSE(xT1, x1_Long)
% RRMSEall_1=[RRMSE_Xi_LL_1(1,1), RRMSE_Xi_LR_1(1,1), RRMSE_Xi_Long_1(1,1)]

%% Evaluate and Plot Results on Testing IC
[tC_LL,xC_LL]=ode45(@(t,x)sparseGalerkin(t,x,Xi_LL,polyorder,usesine),t_span,x05,options); 
[tC_LL,xC_LR]=ode45(@(t,x)sparseGalerkin(t,x,Xi_LR,polyorder,usesine),t_span,x05,options); 
[t1,x1_5]=ode45(@(t,x)sparseGalerkin(t,x,Xi1_Long,polyorder,usesine),t_span,x05,options); 

 figure
hold on
plot(tT5,xT5(:,1),'k','LineWidth',1.95)
plot(tC_LL,xC_LL(:,1),'r--','LineWidth',1.75)
plot(tC_LR,xC_LR(:,1),'g--','LineWidth',1.5)
plot(t1,x1_5(:,1),'b--','LineWidth',1.25)
if eq==1
title('Two-Attractor X-Time Series of True and SINDys on IC5 (Testing Validation)')
legend('True','LL', 'LR','One Long')
end
if eq==2
title('Duffing X-Time Series of True and SINDys on IC5 (Testing Validation)')
legend('True','LL', 'LR','One Long')
end
if eq==3
title('Lorenz X-Time Series of True and SINDys on IC5 (Testing Validation)')
legend('True','II', 'IO','I')
end


ylabel('x')
xlabel('t')

RRMSE_Xi_LL_5=RRMSE(xT5, xC_LL)
RRMSE_Xi_LR_5=RRMSE(xT5, xC_LR)
RRMSE_Xi_Long_5=RRMSE(xT5, x1_5)
% RRMSEall_5=[RRMSE_Xi_LL_5(1,1),RRMSE_Xi_LR_5(1,1), RRMSE_Xi_Long_5(1,1)]

%% Segment Lengths
if n==2
Seglen1=sqrt(diff(xT1(:,1),1).^2+diff(xT1(:,2),1).^2);
arclen1= sum(Seglen1);

    Seglen2=sqrt(diff(xT2(:,1),1).^2+diff(xT2(:,2),1).^2);
arclen2= sum(Seglen2);

    Seglen3=sqrt(diff(xT3(:,1),1).^2+diff(xT3(:,2),1).^2);
arclen3= sum(Seglen3);


 SeglenL=sqrt(diff(xT1_Long(:,1),1).^2+diff(xT1_Long(:,2),1).^2);
arclenL= sum(SeglenL);
else
Seglen1=sqrt(diff(xT1(:,1),1).^2+diff(xT1(:,2),1).^2+diff(xT1(:,3),1).^2);
arclen1= sum(Seglen1);

    Seglen2=sqrt(diff(xT2(:,1),1).^2+diff(xT2(:,2),1).^2+diff(xT2(:,3),1).^2);
arclen2= sum(Seglen2);

    Seglen3=sqrt(diff(xT3(:,1),1).^2+diff(xT3(:,2),1).^2+diff(xT3(:,3),1).^2);
arclen3= sum(Seglen3);


 SeglenL=sqrt(diff(xT1_Long(:,1),1).^2+diff(xT1_Long(:,2),1).^2+diff(xT1_Long(:,3),1).^2);
arclenL= sum(SeglenL);
end


SeglenLimL=Seglen2(end)
SeglenLimR=Seglen3(end)
SeglenLimLong=SeglenL(end)
ArclenLL=arclen1+arclen2
ArclenLR=arclen1+arclen3
arclenL



%%
function Error=RRMSE(xT, xP)
Error=rmse(xT,xP)./std(xT);
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


% function [arclen,seglen] = arclength(X,n)
% 
% if n==2
%  seglen=sqrt(diff(X(:,1),1).^2+diff(X(:,2),1).^2);
% arclen = sum(seglen);
% else
% seglen=sqrt(diff(X(:,1),1).^2+diff(X(:,2),1).^2+diff(X(:,3),1).^2);
% arclen = sum(seglen);
% end 
% end

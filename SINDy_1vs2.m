%% Initialization and System Selection
clear all; clc; close all;
set(0, 'DefaultAxesFontSize', 14);
set(0, 'DefaultTextFontSize', 14);

% Select the system type: 1 = Two-Attractor, 2 = Duffing, 3 = Lorenz
eq = 2;



%% Define System Parameters
    % Lorenz's parameters (chaotic)
        sigma = 10;  
        beta = 8/3;
        rho = 28;

    %Duffing Paraters
        lambda_D=1.3;
        gamma=.5;


% Time Configuration
dt = 0.001;  % Time step
t_final = 25;  % Total time
t_span = 0.001:dt:t_final;  % Time vector

% Initial conditions
    %Two Attractor
        % x01=[-2;-2];  
        % x02=[-1/2;-5];
        %   x03=[1/2;2];
          % x02=[3;-4];

          %Test case the breaks 
% x01=[3;-4];


       % x05=[5;1/25];
       % x05=[2;2];
        
     %Duffing

     
%IC 1 and 2 on different lobs
     % x01=[-1;7/4];  
        % x02=[-3/2;-5];
           x02=[0;-11/3];
          % x04=[3;-4];
%2IC with 1 lobe
     % x01=[3;-4];
      % x02=[-3/2;-5];

%Tight on left
% x01=[-2.3;-0];

%Tight on Right
 x01=[0.79;.75];

       % x05=[-5;2/5];  
     %Lorenz
        % x01=[-2;-2;0];  %right lobe

        %  x01=[2;-5;7];
        % x02=[1/2;2;-2];  %left lobe
        % x03=[-1/2;-5;3];
        % x04=[3;-4;1];
        % 
        % x05=[2;5;15];


%% Integrate System with ode45        
n=2; % number of observables
options = odeset('RelTol',1e-12,'AbsTol',1e-12*ones(1,n));
[tT1,xT1]=ode45(@(t,x)System(t,x, lambda_D, gamma, sigma, beta, rho,eq),t_span,x01,options);
[tT2,xT2]=ode45(@(t,x)System(t,x, lambda_D, gamma, sigma, beta, rho,eq),t_span,x02,options);
 % [tT,xT3]=ode45(@(t,x)System(t,x, lambda_D, gamma, sigma, beta, rho,eq),t_span,x03,options); 
 % [tT,xT4]=ode45(@(t,x)System(t,x, lambda_D, gamma, sigma, beta, rho,eq),t_span,x04,options); 

%% Plot Trajectories
figure
hold on
title('Trajectories of Duffing')
plot(xT1(:,1),xT1(:,2))
plot(xT2(:,1),xT2(:,2))
 % plot(xT3(:,1),xT3(:,2))
% plot(xT4(:,1),xT4(:,2))
grid on
legend('x01','x02','x03','IC4')
xlabel('x')
ylabel('y')
xlim([-5/2,5/2])
ylim([-5/2,5/2])



%% Add Noise to Trajectories

rng(7401) % Seed for reproducibility
%Noise method 1 "+- <= percentage of noise.
        noise_level = 0.005;  % Noise level (0.5%)

       % Generate random noise and apply to both trajectories
        noise1 = (rand(size(xT1)) - 0.5) * 2 * noise_level .* xT1;
        noise2 = (rand(size(xT2)) - 0.5) * 2 * noise_level .* xT2;
        % noise3 = (rand(size(xT3)) - 0.5) * 2 * noise_level .* xT3;
        % noise4 = (rand(size(xT4)) - 0.5) * 2 * noise_level .* xT4;
        

        % Add the noise to the original matrix
        xT1N = xT1 + noise1;
        xT2N=xT2+noise2;
        % xT3N=xT3+noise3;
        % xT4N=xT4+noise4;
    
        % Plot clean vs noisy data
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
%% Define SINDy Library Parameters
polyorder = 3;  % Polynomial library order
usesine = 0; % Disable sine terms

%dimension set
n=n; %preset already for ODE45


%% Define True System Coefficients
 XiT=SystemXiT( lambda_D, gamma, sigma, beta, rho,eq);
poolDataLIST({'x','y'},XiT,n,polyorder,usesine);
% poolDataLIST({'x','y','z'},XiT,n,polyorder,usesine);


%% Compute Savitzky-Golay Derivatives
% Set Savitzky-Golay filter parameters
windowSize = 25;  % Odd number, e.g., 11 points in the window
polynomialOrder = 3;  % Polynomial order, e.g., cubic; dicates how big of a derivate you can take. this can take up to 3rd derivative

% Smooth the signal
% smoothedSignal = sgolayfilt(xT1N, polynomialOrder, windowSize);
% plot(tT,smoothedSignal(:,1),'LineWidth',1.5)
% xlim([6,8])
% Savitzky-Golay differentiation coefficients
[b, g] = sgolay(polynomialOrder, windowSize);

% Calculate the first derivative using the Savitzky-Golay filter
dxT1N=zeros(size(xT1N));
for i=1:n
    dxT1N(:,i)=conv(xT1N(:,i), factorial(1)*g(:,2)/(-dt)^1, 'same');
end
dxT1N=dxT1N((windowSize+1)/2:end-((windowSize-1)/2),:);

dxT2N=zeros(size(xT1N));
for i=1:n
    dxT2N(:,i)=conv(xT2N(:,i), factorial(1)*g(:,2)/(-dt)^1, 'same');
end
dxT2N=dxT2N((windowSize+1)/2:end-((windowSize-1)/2),:);

% dxT3N=zeros(size(xT1N));
% for i=1:n
%     dxT3N(:,i)=conv(xT3N(:,i), factorial(1)*g(:,2)/(-dt)^1, 'same');
% end
% dxT3N=dxT3N((windowSize+1)/2:end-((windowSize-1)/2),:);
% 
% dxT4N=zeros(size(xT1N));
% for i=1:n
%     dxT4N(:,i)=conv(xT4N(:,i), factorial(1)*g(:,2)/(-dt)^1, 'same');
% end
% dxT4N=dxT4N((windowSize+1)/2:end-((windowSize-1)/2),:);

%% Align and Combine Derivative and State Data
 
 xT1N=xT1N((windowSize+1)/2:end-((windowSize-1)/2),:);
 xT2N=xT2N((windowSize+1)/2:end-((windowSize-1)/2),:);
 % xT3N=xT3N((windowSize+1)/2:end-((windowSize-1)/2),:);
 % xT4N=xT4N((windowSize+1)/2:end-((windowSize-1)/2),:);

 xTC=[xT1N; xT2N];
   % xTC=[xT1N; xT2N; xT3N;xT4N];

dxTC=[dxT1N; dxT2N];
% dxTC=[dxT1N; dxT2N;dxT3N;dxT4N];


%% Build Library of Nonlinear Terms
Theta_C= poolData(xTC,n,polyorder,usesine);
Theta_1= poolData(xT1N,n,polyorder,usesine);


%% Perform Sparse Regression
lambda = 0.25;      % lambda is our sparsification threshold.

iters=25;
XiC = sparsifyDynamics(Theta_C,dxTC,lambda,n,iters);
poolDataLIST({'x','y'},XiC,n,polyorder,usesine);
% poolDataLIST({'x','y','z'},XiC,n,polyorder,usesine);
Xi1 = sparsifyDynamics(Theta_1,dxT1N,lambda,n,iters);
poolDataLIST({'x','y'},Xi1,n,polyorder,usesine);
% poolDataLIST({'x','y','z'},Xi1,n,polyorder,usesine);

%% Results
%% Evaluate on Training Trajectory (x01)
[tC,xC_1]=ode45(@(t,x)sparseGalerkin(t,x,XiC,polyorder,usesine),t_span,x01,options); 
 [t1,x1_1]=ode45(@(t,x)sparseGalerkin(t,x,Xi1,polyorder,usesine),t_span,x01,options); 

Sparsity_XiC=Sparsity(XiT,XiC)
Sparsity_Xi1=Sparsity(XiT,Xi1)
Sparsity_all=[Sparsity_XiC, Sparsity_Xi1]

figure
hold on
title('Lorenz X-Time Series of True and SINDys on IC1 (Training Validation)')
plot(tT1,xT1(:,1),'k','LineWidth',1.75)
plot(tC,xC_1(:,1),'r--','LineWidth',1.5)
plot(t1,x1_1(:,1),'b--','LineWidth',1.25)

legend('True','Multi-Trajectory','Single Trajectory')
ylabel('x')
xlabel('t')

RRMSE_XiC_1=RRMSE(xT1, xC_1);
RRMSE_Xi1_1=RRMSE(xT1, x1_1);
RRMSEall_1=[RRMSE_XiC_1, RRMSE_Xi1_1]

%% Evaluate on Test Trajectory (x05)
[tT5,xT5]=ode45(@(t,x)System(t,x, lambda_D, gamma, sigma, beta, rho,eq),t_span,x05,options);
[tC,xC_5]=ode45(@(t,x)sparseGalerkin(t,x,XiC,polyorder,usesine),t_span,x05,options); 
[t1,x1_5]=ode45(@(t,x)sparseGalerkin(t,x,Xi1,polyorder,usesine),t_span,x05,options); 

 figure
hold on
title('Lorenz X-Time Series of True and SINDys on IC5 (Testing Validation)')
plot(tT5,xT5(:,1),'k','LineWidth',1.75)
plot(tC,xC_5(:,1),'r--','LineWidth',1.5)
plot(t1,x1_5(:,1),'b--','LineWidth',1.25)

legend('True','Multi-Trajectory','Single Trajectory')
ylabel('x')
xlabel('t')

RRMSE_XiC_5=RRMSE(xT5, xC_5)
RRMSE_Xi1_5=RRMSE(xT5, x1_5);
RRMSEall_5=[RRMSE_XiC_5, RRMSE_Xi1_5]

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
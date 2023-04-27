clear all

clf

% Input Parameters
a=0.005;
b=0.004;
Dc=1e-6; % [m]
Vl=1e-4; % Loading Rate [m/s]
Vini=1e-5; % Initial Velocity [m/s]
ThetaI=Dc/Vini; % Initial State - set Start from steady state
NormalStress=1e6; % 2Mpa
V0=1e-9; % Reference velocity [m/s]
Mass=1000; % per unit area [kg/m^2] (for example, 3000=30kg on 0.1 m by 0.1 m]
K=1e9;%0.8*Kc_QS; % Set stiffness 0.8 of critical stiffness
Friction0=0.6; % reference friction

FrictionI=Friction0+a*log(Vini/V0)+b*log(ThetaI*V0/Dc); % Initial friction
Xl_Ini=FrictionI*NormalStress/K; % initial load point


% Time Step Control KJ
Dt=1e-4; % Time step [second]
TotalTime=2; % Total time [second]
TotalStep=round(TotalTime/Dt); % Total steps

% RK4 TimeStep
dt = 1e-6; % time step ODE
tspan = [0, TotalTime]; % time span ODE Solver



% Convergence control KJ Method
V_eps=1e-7 % Convergence criterion in NR
% DV=1e-15   % Denominator of NR - this will be changed with tested Velocity
Theta_eps=1e-6; % Convergence criterion for Theta update (implicit only)
DelTheta=1e-5; % Denominator of NR for Theta update (implicit only)


%%%%%%%%%%%%%%%%%%% KJ's Method %%%%%%%%%%%%%%%%%%%%%
tic
% Simulation begins
XlOld=Xl_Ini;
Omega=sqrt(K/Mass);
Step=0;
for i=1:TotalStep
    
    Xl=XlOld+Dt*Vl; % Load point displacement

    
    if i==1
        Friction=FrictionI;
        DispOld=0;
        Disp=0;
        Theta=ThetaI;
        ThetaOld=Theta;
        VOld=Vini;
    end
    
    % Newton Rhapson Begins
    VDiff=10; % Arbitrary for initiation
    Iteration=0; % Number of iteration
    BREAK=0; % This only used when convergence is hard to made
%     V=1e-11; % Arbitrary initial velocity (good to be small to pick up small V)
    V=VOld/1e10; % Arbitrary initial velocity (good to be small to pick up small V)
    while abs(VDiff-1)>V_eps
        Iteration=Iteration+1;
        VTest=V; % Velocity tested in this NR iteration
        DV=V/1000; % Denominator of NR - changes with tested velocity
        
        % Finding Initial Value
        VOldIter=V;

%         Theta=(ThetaOld-Dc/V)*exp(-V/Dc*Dt)+Dc/V;
        Theta=(ThetaOld+Dt)/(1+V*Dt/Dc); % Deterich Evolution        
        Friction=Friction0+b*log(V0*Theta/Dc)+a*log(V/V0) ; % Equation (15) Rate and State Friction
        F=Xl-Friction*NormalStress/K;
        Disp=(DispOld-F)*cos(Omega*Dt)+(VOld/Omega)*sin(Omega*Dt)+F; % Equation (9)
        V=(Disp-DispOld)/Dt*2-VOld; % Equation (12)
        
        FOriginal=VOldIter-V; % We are testing this Newton Rhapson Function. Lets send this to zero
        
        % Finding deviated value for NR
        V=VOldIter+DV;
        VOldIter=V;
%         Theta=(ThetaOld-Dc/V)*exp(-V/Dc*Dt)+Dc/V;        
        Theta=(ThetaOld+Dt)/(1+V*Dt/Dc); % Deterich Evolution                
        Friction=Friction0+b*log(V0*Theta/Dc)+a*log(V/V0) ; % Equation (15) Rate and State Friction
        F=Xl-Friction*NormalStress/K;
        Disp=(DispOld-F)*cos(Omega*Dt)+(VOld/Omega)*sin(Omega*Dt)+F; % Equation (9)
        V=(Disp-DispOld)/Dt*2-VOld; % Equation (12)
        
        NRF=VOldIter-V; % Recalculate NR testing function with this velocity
        
        DF=(NRF-FOriginal)/DV; % tangent of the NR function
        V=VTest-FOriginal/DF; % Update velocity
        
        if V<0; BREAK=1; % Only used when convergence is failed
            break; end
        
        VDiff=abs(VTest/V); % Calculate the Convergence criterion
        if Iteration>100; V_eps=V_eps*2 % Just in case it is too hard to be converged
        end
    end % End of NR iteration
    
    
    if BREAK==1  % only if we could not get convergence due to very small velocity
        V=0; % set it just zero
        Theta=ThetaOld+Dt;
        Disp=DispOld;
    end
    
    
    if rem(i,1)==0 % save the data in every 1 steps
        T=i*Dt;
        %[T,log10(V)]
        Step=Step+1;
        VHistory(Step)=V; % velocity
        DispHistory(Step)=Disp; % Dispolacement
        ThetaHistory(Step)=Theta; % State variable
        Time(Step)=i*Dt; % Time
    end
    
    DispOld=Disp;
    ThetaOld=Theta;
    VOld=V;
    XlOld=Xl;
end

fprintf("KJ's solver ")
toc

% Plots
figure(1)
hold on
set(gcf, 'color', 'w')
set(gca,'fontsize', 13)
ylabel('Velocity')
xlabel('Time (s)')
plot(Time,VHistory, 'k', 'LineWidth',2)
set(gca, 'YScale', 'log')
box on
drawnow









% ode = @(t, y) [ 1/(NormalStress*a/y(1) + Eta) * ( K*(Vl-y(1)) - NormalStress*b/y(2)*(1-y(1)*y(2)/Dc)); 1-y(1)*y(2)/Dc];

ode = @(t, y) [y(2); K/Mass*(t*Vl + Xl_Ini - y(1)) - (Friction0 + a*log(y(2)/V0) + b*log(V0*y(3)/Dc))*NormalStress/Mass; 1-y(2)*y(3)/Dc];


% Initialize the solution
t = tspan(1):dt:tspan(2);
N = length(t);
Disp_RK = zeros(N, 1);
Vel_RK = zeros(N, 1);
Theta_RK = zeros(N, 1);

% Set the initial conditions
Disp_RK(1) = 0;
Vel_RK(1) = Vini;
Theta_RK(1) = ThetaI;

tic
for i = 1:N-1
    % Calculate the derivatives at the current time step
    k1 = ode(t(i), [Disp_RK(i); Vel_RK(i); Theta_RK(i)]);
    k2 = ode(t(i) + dt/2, [Disp_RK(i) + dt*k1(1)/2; Vel_RK(i) + dt*k1(2)/2; Theta_RK(i) + dt*k1(3)/2]);
    k3 = ode(t(i) + dt/2, [Disp_RK(i) + dt*k2(1)/2; Vel_RK(i) + dt*k2(2)/2; Theta_RK(i) + dt*k2(3)/2]);
    k4 = ode(t(i) + dt, [Disp_RK(i) + dt*k3(1)/2; Vel_RK(i) + dt*k3(2); Theta_RK(i) + dt*k3(3)]);
    
    % Update the solution at the next time step
    Disp_RK(i+1) = Disp_RK(i) + dt*(k1(1) + 2*k2(1) + 2*k3(1) + k4(1))/6;
    Vel_RK(i+1) = Vel_RK(i) + dt*(k1(2) + 2*k2(2) + 2*k3(2) + k4(2))/6;
    Theta_RK(i+1) = Theta_RK(i) + dt*(k1(3) + 2*k2(3) + 2*k3(3) + k4(3))/6;
end
fprintf("RK4 ")
toc

% Plot the displacement of the damped oscillator over time
figure(1)
% Plot the displacement of the damped oscillator over time
plot(t, Vel_RK,'r--', 'LineWidth',2)
xlabel('Time')
ylabel('V')
set(gca, 'YScale', 'log')


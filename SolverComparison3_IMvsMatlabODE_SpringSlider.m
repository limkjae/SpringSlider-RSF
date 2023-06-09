clear all

clf

% Input Parameters
a=0.005;
b=0.007;
Dc=1e-6; % [m]
Vl=1e-4; % Loading Rate [m/s]
Vini=1e-6; % Initial Velocity [m/s]
ThetaI=Dc/Vini; % Initial State - set Start from steady state
NormalStress=1e6; % 
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

% Matlab ODE Solver MaxTimeStep
Maxdt_ODE = 1e-4; % time step ODE
tspan = [0, TotalTime]; % time span ODE Solver



% Convergence control KJ Method
V_eps=1e-7 % Convergence criterion in NR
% DV=1e-15   % Denominator of NR - this will be changed with tested Velocity
Theta_eps=1e-6; % Convergence criterion for Theta update (implicit only)
DelTheta=1e-5; % Denominator of NR for Theta update (implicit only)

DorR=1; % 1:Deterich, 2:Ruina
EorI=2; % Theta update 1: explicit, 2: implicit



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
        Friction=Friction0+b*log(V0*Theta/Dc)+a*log(V/V0); % Equation (15) Rate and State Friction
        F=Xl-Friction*NormalStress/K;
        Disp=(DispOld-F)*cos(Omega*Dt)+(VOld/Omega)*sin(Omega*Dt)+F; % Equation (9)
        V=(Disp-DispOld)/Dt*2-VOld; % Equation (12)
        
        FOriginal=VOldIter-V; % We are testing this Newton Rhapson Function. Lets send this to zero
        
        % Finding deviated value for NR
        V=VOldIter+DV;
        VOldIter=V;
%         Theta=(ThetaOld-Dc/V)*exp(-V/Dc*Dt)+Dc/V;        
        Theta=(ThetaOld+Dt)/(1+V*Dt/Dc); % Deterich Evolution                
        Friction=Friction0+b*log(V0*Theta/Dc)+a*log(V/V0); % Equation (15) Rate and State Friction
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
        Time(Step)=T; % Time
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
% plot(Time,RASFricHistory,'r')
box on
drawnow
DT_History_IM=Time(2:end)-Time(1:end-1)
% plot(DT_History_IM)





%%%%%%%%%%%%%%%% Matlab ODE %%%%%%%%%%%%%%%%%%%%%
tic
ode = @(t, y) [y(2); K/Mass*(t*Vl + Xl_Ini - y(1)) - (Friction0 + a*log(y(2)/V0) + b*log(V0*y(3)/Dc))*NormalStress/Mass; 1-y(2)*y(3)/Dc];
[t, y] = ode23s(ode, tspan, [0, Vini, ThetaI], odeset('MaxStep', Maxdt_ODE));
fprintf("Matlab ODE ")
toc
if isreal(y)==0
   fprintf("Matlab ODE Failed \n") 
end



% ode = @(t, y) [y(2); K/Mass*(t*Vl + Xl - y(1)) - (Friction0 + (a-b)*log(y(2)/V0))*NormalStress/Mass];
% [t, y] = ode45(ode, tspan, [Vini, 0], odeset('MaxStep', dt));

figure(1)
% Plot the displacement of the damped oscillator over time
plot(t, y(:,2),'r--', 'LineWidth',2)
% xlabel('Time')
% ylabel('V')
set(gca, 'YScale', 'log')

% 
% 
% DT_History_Matlab=t(2:end)-t(1:end-1)
% figure(2)
% plot(DT_History_IM,'k')
% hold on
% plot(DT_History_Matlab,'r--')

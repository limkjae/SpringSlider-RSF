%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Spring Slider with rate and state friction solver Matlab version.\
%%% Written by Kyungjae (KJ) Im
%%% For details - 
%%% Im, K., Elsworth, D., Marone, C., & Leeman, J. (2017). The impact of 
%%% frictional healing on stick-slip recurrence interval and stress drop: 
%%% Implications for earthquake scaling. Journal of Geophysical Research: 
%%% Solid Earth. https://doi.org/10.1002/2017JB014476 
%%% (Section 3. Numerical method of Stick-slip simulation)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all

clf

% Input Parameters
a=0.003;
b=0.008;
Dc=1e-6; % [m]
Vl=3e-5; % Loading Rate [m/s]
Vini=3e-7; % Initial Velocity [m/s]
ThetaI=Dc/Vini; % Initial State - set Start from steady state
NormalStress=40e6; % 2Mpa
V0=3e-6; % Reference velocity [m/s]
Mass=2000; % per unit area [kg/m^2] (for example, 3000=30kg on 0.1 m by 0.1 m]
K=90.5e9;%0.8*Kc_QS; % Set stiffness 0.8 of critical stiffness
Friction0=0.7; % reference friction

FrictionI=Friction0+a*log(Vini/V0)+b*log(ThetaI*V0/Dc); % Initial friction
Xl_Ini=FrictionI*NormalStress/K; % initial load point


% Time Step Control KJ
Dt = 1e-05 ; % Time step [second]
if Dt>sqrt(Mass/K)/10
    fprintf("!!! Warning!!!!!!!!!!!! The time step may be too large. \n!!! We recommend Dt < sqrt(Mass/K)/10 \n")
    fprintf("!!! Consider Dt < %s \n", sqrt(Mass/K)/10)
end
TotalTime=10; % Total time [second]
TotalStep=round(TotalTime/Dt); % Total steps


% Convergence control KJ Method
V_eps=1e-7; % Convergence criterion in NR
% DV=1e-15   % Denominator of NR - this will be changed with tested Velocity
Theta_eps=1e-6; % Convergence criterion for Theta update (implicit only)
DelTheta=1e-5; % Denominator of NR for Theta update (implicit only)


%%%%%%%%%%%%%%%%%%% KJ's Method %%%%%%%%%%%%%%%%%%%%%
% tic
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
        ShearStressHistory(Step)=Friction*NormalStress;
        FrictionHistory(Step)=Friction;
        
        Time(Step)=i*Dt; % Time
    end
    
    DispOld=Disp;
    ThetaOld=Theta;
    VOld=V;
    XlOld=Xl;
end

% fprintf("KJ's solver ")

% toc
% 
% if Dt>sqrt(Mass/K)/10
%     fprintf("!!! Consider Dt < %s \n", sqrt(Mass/K)/10)
% end

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

% 
% % Plots
% figure(2)
% hold on
% set(gcf, 'color', 'w')
% set(gca,'fontsize', 13)
% % ylabel('Velocity')
% xlabel('Time (s)')
% plot(Time,ShearStressHistory, 'k', 'LineWidth',2)
% % set(gca, 'YScale', 'log')
% box on
% drawnow




%Model name: model of Memristor based on Theoretical formulas
%           this model was written by Dmitry Fliter and Keren Talisveyberg 
%           Technion Israel institute of technology EE faculty December 2011


model = 0;   % define the model 0 - Linear Ion Drift ; 1 - Simmons Tunnel Barrier; 2 - Team model ; 3 - Nonlinear Ion Drift model ; 4 - VTEAM
win   = 3;   % define the window type :  0 - No window; 1 - Jogelkar window ; 2 - Biolek window ; 3 - Prodromakis window ; 4- Kvatinsky window (Team model only) 
iv    = 0;   % IV_relation=0 means linear V=IR. IV_relation=1 means nonlinear V=I*exp{..}  

%Genaral parameters
num_of_cycles = 10;
amp = 1;
freq = 0.3;
w_init = 0.5; % the initial state condition [0:1] 
D = 11;
V_t = 0;
P_coeff = 3;
J = 1;
Roff = 16e3;
Ron = 130;

%Linear Ion Drift parameters 
uV=1e-15;                                %%dopants mobility

%Nonlinear Ion Drift parametrs
beta = 9;
a = 4;
c = 0.01;
n = 14;
q = 13;
g = 4;
alpha = 7;

%% section 1: %% Linear Ion Drift model.

tspan=[0 num_of_cycles/freq];                       %%time length of the simulation
points=2e5;                              %%number of sampling points
W0=w_init*D;                            %define the initial value of W
tspan_vector = linspace(tspan(1),tspan(2),points);         % Create vector of initial values
V = amp*sin(freq*2*pi*tspan_vector);
current=zeros(size(tspan_vector));
W=zeros(size((tspan_vector)));
W_dot=zeros(size((tspan_vector)));
delta_t=tspan_vector(2)-tspan_vector(1);                        %%define the step size

W(1)=W0;                                                 %% initiliaze the first W vetor elemnt to W0 - the initial condition
for i=2:length(tspan_vector)
         
  % case this is Prodromakis window
  W_dot(i)=a*(V(i)^q);
  W(i)=W(i-1)+W_dot(i)*delta_t*(J*(1-((W(i-1)/D-0.5)^2+0.75)^P_coeff));

  % correct the w vector according to bounds [0 D]
    if W(i) < 0
        W(i) = 0;
        W_dot(i)=0;
    elseif W(i) > D
        W(i) = D;
        W_dot(i)=0;
    end

current(i)=W(i)^n*beta*sinh(alpha*V(i))+c*(exp(g*V(i))-1);
end

R=Ron*W/D+Roff*(1-W/D); %%this parameter might be useful for debug
figure(1);
plot(V(20e3:end),current(20e3:end));
hold on;
title('I-V curve');
xlabel('V[volt]');
ylabel('I[amp]');
plot([0 0], [-6 6], 'black');  %x-axis
plot([-1.2 1.2], [0 0], 'black');  %y-axis

figure(2);
plot(tspan_vector,W/D,'r');
title('W/D as func of time');
xlabel('time[sec]');
legend('W/D');

figure(3);
plot(tspan_vector,R,'black');
title('R as func of time');
xlabel('time[sec]');
legend('R');

%%  Nonlinear Ion Drift model section C with capacitor

%variable overriding
amp = 1;
freq = 1.3;


points=5e5;                              %%number of sampling points
tspan=[0 num_of_cycles/freq];
tspan_vector = linspace(tspan(1),tspan(2),points);         % Create vector of initial values
W0=w_init*D;                            %define the initial value of W
t = linspace(tspan(1),tspan(2),points);
delta_t = t(2) - t(1);
V = amp*square(2*pi*freq*t);
Vm = zeros(size(t));
cap = 0.01;
Vc = zeros(size(t));
Vc_max = amp;
Vc_min = -amp;
Vc_max_zero = 0;
Vc_min_zero = 0;
G_effective = 0;
R_effective_vec = zeros(size(t));
W=zeros(size((t)));
W_dot=zeros(size((t)));
current=zeros(size((t)));

W(1) = W0;
for i=2:points
    W_dot(i)=a*(V(i)^q);
    W(i)=W(i-1)+W_dot(i)*delta_t*(J*(1-((W(i-1)/D-0.5)^2+0.75)^P_coeff));
    

    G_effective = (W(i-1)^n * beta * alpha * cosh(alpha * V(i-1))) + (c * g * exp(g * V(i-1))); %remove
    R_effective = 1 / G_effective;
    R_effective_vec(i) = R_effective;
    
    if (V(i) > 0)
        Vc(i) = Vc_max - (Vc_max - Vc_min_zero) * exp(-mod(t(i), 1/(2*freq)) / (R_effective*cap));
        if Vc_max_zero < Vc(i)
            Vc_max_zero = Vc(i);
        end
    elseif (V(i) < 0)
        Vc(i) = Vc_min - (Vc_min - Vc_max_zero) * exp(-mod(t(i), 1/(2*freq)) / (R_effective*cap)); %remove
        if Vc_min_zero > Vc(i) %remove
            Vc_min_zero = Vc(i); % remove
        end
    end
    Vm(i)=V(i)-Vc(i); 

    % correct the w vector according to bounds [0 D]
    if W(i) < 0
        W(i) = 0;
        W_dot(i)=0;
    elseif W(i) > 1
        W(i) = 1;
        W_dot(i)=0;
    end
    current(i)=W(i)^n*beta*sinh(alpha*Vm(i))+c*(exp(g*Vm(i))-1); % remove
end 

figure(1);
plot(Vm(20e2:end),current(20e2:end));
title('I-V curve');
xlabel('V[volt]');
ylabel('I[amp]');

figure(2);
plot(t,W,'r');
title('W/D');
xlabel('time[sec]');
legend('W/D');

figure(3);
plot(t,Vm,'r');
title('Vm(t)');
xlabel('time[sec]');
legend('Vm');

figure(4);
plot(t,V);
title('V(t)');
xlabel('time[sec]');
legend('V');

figure(5);
plot(t,R_effective_vec);
title('M(t)');
xlabel('time[sec]');
legend('M');

figure(6);
plot(t,current);
title('I(t)');
xlabel('time[sec]');
legend('I');

disp(Vc_min_zero);
disp(Vc_max_zero);



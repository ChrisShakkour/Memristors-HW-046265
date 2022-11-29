%% section b
Qb = 1;
Rb = 1;
t = linspace(0,20, 1000);
Ib = sin(t);
qb = cos(t) + Qb;
Mb = Rb./sqrt(1+((qb.*qb)/Qb));
phib = Qb*Rb*log((qb/Qb)+sqrt(1+((qb.*qb)/Qb^2)));
Vb = Mb .* Ib;

figure(1);
plot(t,Ib);
hold on;
plot(t,Vb);
title('V(t),I(t) plot');
xlabel('T[sec]');
ylabel('V[v], I[amp]');
legend('I(t)', 'V(t)');

figure(2);
plot(Vb,Ib);
title('V-I curve plot');
xlabel('V[v]');
ylabel('I[amp]');

figure(3);
plot(t,Mb);
title('Memristance plot');
xlabel('T[sec]');
ylabel('M(t)');


%% section c
Qc = 1;
Rc = 1;
K = 1 + sqrt(2);
Vc = sin(t);
qc = (K^2 * exp(2*cos(t)) - 1) ./ (2 * K * exp(cos(t)));
Mc = Rc./sqrt(1+((qc.*qc)/Qc));
phic = Qc*Rc*log((qc/Qc)+sqrt(1+((qc.*qc)/Qc^2)));
%V = diff(phi);
Ic = Vc ./ Mc;

figure(1);
plot(t, Vc, t, Ic);
title('V(t)/I(t)');
legend('v[V]','I[A]');
xlabel('time[sec]');
ylim([-2 2]);

figure(2);
plot(Vc, Ic);
title('I-V curve plot');
xlabel('v[V]');
ylabel('I[A]');
ylim([-2 2]);

figure(3);
plot(t,Mc);
title('Memristance plot');
xlabel('T[sec]');
ylabel('M(t)');

%% section e

w = 0.3; %pulse width
d = 0:w*5:20; %delay vector
Ie = pulstran(t,d,'rectpuls',w);
Ie_p = Ie * 0.1; %train of positive pulses
Ie_n = Ie_p * (-1);%train of negative pulses
qe_p = zeros(1, length(Ie_p));
qe_n = zeros(1, length(Ie_n));
for i = 1:length(Ie_p)
    qe_p(i) = sum(Ie_p(1:i));
    qe_n(i) = sum(Ie_n(1:i));
end
Me_n = Rc./sqrt(1+((qe_n.*qe_n)/Qc));
Me_p = Rc./sqrt(1+((qe_p.*qe_p)/Qc));

figure(2);
subplot(2,2,1);
plot(t, qe_n, t, qe_p);
legend('q(t) for negative pulses','q(t) for positive pulses');
xlabel('time[sec]');
ylabel('Q[c]');
title('q(t)');

subplot(2,2,2);
plot(t,Me_p);
ylabel('M[Ohm]');
title('Memristance(t)');

subplot(2,2,3);
plot(t, Ie_p);
title('I(t) - positive pulses');
xlabel('time[sec]');
ylabel('I[A]');

subplot(2,2,4);
plot(t, Ie_n);
title('I(t) - negative pulses');
xlabel('time[sec]');
ylabel('I[A]');





clear all; clc;close all;

syms t I3 I2 I1

g=9.8;
m3 = 50;
m2 = 50;
m1 = 50;
d1 = 50;

lc1 = 0.5;
a2 =1 ;
lc2 =0.5;
a3 = 1;
lc3 = 0.5;

I3 = 10;
I2 = 10;
I1 = 10;


q1(t) = sin(t);
q2(t) = cos(2*t);
q3(t) = sin(3*t);

q1_dot(t) = cos(t);
q2_dot(t) = -2*sin(2*t);
q3_dot(t) = 3*cos(3*t);

d_2dot(t) = 0*t;
q2_2dot(t) = 0*t;
q3_2dot(t) = 0*t;


R1_2 = rotz(q2(t));
R2_3 = rotz(q3(t));

r1_0 = [0;d1;0];
r1_1 = [0;lc1;0];
r2_1 = [a2;0;0];
r2_2 = [lc2;0;0];
r3_2 = [a3;0;0];
r3_3 = [lc3;0;0];

w1 = [0;0;0];; %omega
w1_dot = 0; %omega dot

w2 = [0;0;q2_dot(t)];
w2_dot = [0;0;q2_2dot(t)];

w3 = [0;0;q2_dot+q3_dot(t)];
w3_dot = [0;0;q2_2dot(t)+q3_2dot(t)];


P_2dot_1 = [-g;d_2dot(t);0];
P_2dot_c1 = P_2dot_1;

P_2dot_2 = R1_2.'*P_2dot_1+cross(w2_dot,r2_1)+cross(w2,cross(w2,r2_1));
P_2dot_c2 = P_2dot_2 +cross(w2_dot,r2_2)+cross(w2,cross(w2,r2_2));

P_2dot_3 = R2_3.'*P_2dot_2+cross(w3_dot,r3_2)+cross(w3,cross(w3,r3_2));
P_2dot_c3 = P_2dot_3 +cross(w3_dot,r3_3)+cross(w3,cross(w3,r3_3));

%backward 
f3_3 = m3*P_2dot_c3;
Tau_3(t) = I3*w3_dot+cross(w3,I3*w3)-cross(f3_3,r3_3);

f2_2 = R2_3*f3_3+m2*P_2dot_c2;
Tau_2(t) = R2_3*Tau_3-cross(f2_2,r2_2)-cross(R2_3*f3_3,r2_1-r2_2)+I2*w2_dot+cross(w2,I2*w2);

f1_1 = R1_2*f2_2+m1*P_2dot_c1;
Tau_1(t) = R1_2*Tau_2-cross(f1_1,r1_1)-cross(R1_2*f2_2,r1_0-r1_1)+I1*w1_dot+cross(w1,I1*w1);

a = Tau_1(t);
b = Tau_2(t);
c = Tau_3(t);


% Set latex interpreter for plots text and legends
set(0, 'DefaultTextInterpreter', 'latex') 
set(0, 'DefaultLegendInterpreter', 'latex')

figure;
subplot(3,1,1);
fplot(a(3), [0 2*pi],'LineWidth',1.5,'color','red')
title('Coriolis of 1st Joint');
xlabel('Time $t$ [s]');
ylabel('Torque $\tau$ [$N.m$]');
grid

subplot(3,1,2);
fplot(b(3), [0 2*pi],'LineWidth',1.5,'color','blue')
title('Coriolis of 2nd Joint');
xlabel('Time $t$ [s]');
ylabel('Torque $\tau$ [$N.m$]');
grid

subplot(3,1,3);
fplot(c(3), [0 2*pi],'LineWidth',1.5,'color','green')
title('Coriolis of 1st Joint');
xlabel('Time $t$ [s]');
ylabel('Torque $\tau$ [$N.m$]');
grid



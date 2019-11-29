clear all; clc;

syms q2 q3 d1 lc1 a2 lc2 a3 lc3 m3 m2 m1 d_2dot q2_dot q2_2dot g q3_dot q3_2dot I3 I2 I1


R1_2 = rotz(q2);
R2_3 = rotz(q3);

r1_0 = [0;d1;0];
r1_1 = [0;lc1;0];
r2_1 = [a2;0;0];
r2_2 = [lc2;0;0];
r3_2 = [a3;0;0];
r3_3 = [lc3;0;0];

w1 = [0;0;0];; %omega
w1_dot = 0; %omega dot

w2 = [0;0;q2_dot];
w2_dot = [0;0;q2_2dot];

w3 = [0;0;q2_dot+q3_dot];
w3_dot = [0;0;q2_2dot+q3_2dot];


P_2dot_1 = [-g;d_2dot;0];
P_2dot_c1 = P_2dot_1;

P_2dot_2 = R1_2.'*P_2dot_1+cross(w2_dot,r2_1)+cross(w2,cross(w2,r2_1));
P_2dot_c2 = P_2dot_2 +cross(w2_dot,r2_2)+cross(w2,cross(w2,r2_2));

P_2dot_3 = R2_3.'*P_2dot_2+cross(w3_dot,r3_2)+cross(w3,cross(w3,r3_2));
P_2dot_c3 = P_2dot_3 +cross(w3_dot,r3_3)+cross(w3,cross(w3,r3_3));

%backward 
f3_3 = m3*P_2dot_c3;
Tau_3 = I3*w3_dot+cross(w3,I3*w3)-cross(f3_3,r3_3)

f2_2 = R2_3*f3_3+m2*P_2dot_c2;
Tau_2 = R2_3*Tau_3-cross(f2_2,r2_2)-cross(R2_3*f3_3,r2_1-r2_2)+I2*w2_dot+cross(w2,I2*w2)

f1_1 = R1_2*f2_2+m1*P_2dot_c1;
Tau_1 = R1_2*Tau_2-cross(f1_1,r1_1)-cross(R1_2*f2_2,r1_0-r1_1)+I1*w1_dot+cross(w1,I1*w1)




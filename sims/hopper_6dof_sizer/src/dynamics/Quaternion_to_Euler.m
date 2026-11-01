function [phi, theta, psi] = Quaternion_to_Euler(q0,q1,q2,q3)



L_11 = (q0^2+q1^2)-q2^2-q3^2;
L_12 = 2*((q1*q2)+q0*q3);
L_13 = 2*((q1*q3)-q0*q2); 

L_21 = 2*((q1*q2)-(q0*q3));
L_22 = (q0^2-q1^2+q2^2-q3^2);
L_23 = 2*((q2*q3)+(q0*q1));

L_31 = 2*((q1*q3)+(q0*q2));
L_32 = 2*((q2*q3)-q0*q1);
L_33 = q0^2-q1^2-q2^2+q3^2;

L_ba = [L_11 L_12 L_13;
        L_21 L_22 L_23;
        L_31 L_32 L_33];

phi = atan2(L_ba(2,3),L_ba(3,3));
theta = -asin(L_ba(1,3));
psi = atan2(L_ba(1,2),L_ba(1,1));

end
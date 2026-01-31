function [A,B]  = HopperLinearization_3DOF(theta,T,delta)


%% JACOBIANS: A and B
m=10;
l=.2;
I = 1;

A = [0 0 1 0 0 0;
     0 0 0 1 0 0;
     0 0 0 0 (-T/m)*sin(theta+delta) 0;
     0 0 0 0 (T/m)*cos(theta+delta) 0;
     0 0 0  0 0 1;
     0 0 0 0 0 0];

B = [0 0;
    0 0;
    (1/m)*cos(theta+delta), (-T/m)*sin(theta+delta);
    (1/m)*sin(theta+delta), (T/m)*cos(theta+delta);
    0 0;
    (l/I)*sin(delta), ((l*T)/I)*cos(delta)];
end



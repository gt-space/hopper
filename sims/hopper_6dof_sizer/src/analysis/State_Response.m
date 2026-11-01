function [X,phi,theta,psi,t,xyz_ts,phi_ts,theta_ts,psi_ts] = State_Response(sys_IC,tresp,disturbed_state,disturbance)

x0 = zeros(13,1);

for i = 1:length(disturbed_state)
x0(disturbed_state(i)) =disturbance(i);
end

[~,t,X] = initial(sys_IC,x0,tresp);

for i = 1:length(t)
[phi(i),theta(i),psi(i)] = Quaternion_to_Euler(X(i,10),X(i,11),X(i,12),X(i,13));
end

phi = rad2deg(phi); theta = rad2deg(theta); psi = rad2deg(psi);
xyz_ts = timeseries(X(:,[1:3]), t);
phi_ts = timeseries(phi', t);
theta_ts = timeseries(theta', t);
psi_ts = timeseries(psi', t);

end
function [innovation, measurementJacobian] = gps_velocity_update_model(nominalState, measurement, params)
%GPS_VELOCITY_UPDATE_MODEL  Independent 3-DOF GPS velocity update.
velocityEstimate = nominalState(4:6);
angularVelocityEstimate = nominalState(7:9);
quaternionEstimate = nominalState(10:13);
r_B = params.r_b(:);
C_IB = quaternion_to_rotation(quaternionEstimate);

predictedVelocity = velocityEstimate + C_IB * cross(angularVelocityEstimate, r_B);
innovation = predictedVelocity - measurement.vel(:);

measurementJacobian = zeros(3,19);
measurementJacobian(:,4:6) = eye(3);
measurementJacobian(:,7:9) = -C_IB * skew(r_B);
measurementJacobian(:,10:12) = -skew(C_IB * cross(angularVelocityEstimate, r_B));
end

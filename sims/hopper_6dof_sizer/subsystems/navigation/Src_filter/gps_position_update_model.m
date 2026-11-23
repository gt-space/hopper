function [innovation, measurementJacobian] = gps_position_update_model(nominalState, measurement, params)
%GPS_POSITION_UPDATE_MODEL  Independent 3-DOF GPS position update.
positionEstimate = nominalState(1:3);
angularVelocityEstimate = nominalState(7:9);
quaternionEstimate = nominalState(10:13);
r_B = params.r_b(:);
C_IB = quaternion_to_rotation(quaternionEstimate);

predictedPosition = positionEstimate + C_IB * r_B;
innovation = predictedPosition - measurement.pos(:);

measurementJacobian = zeros(3,19);
measurementJacobian(:,1:3) = eye(3);
measurementJacobian(:,10:12) = -skew(C_IB * r_B);
end

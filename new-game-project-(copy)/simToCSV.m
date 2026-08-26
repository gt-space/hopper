simmulation = sim("hopper_6dof_NED_v2.slx")
csvArray = array2table([simmulation.tout,simmulation.x.data,simmulation.y.data,simmulation.z.data,simmulation.thrust.data,simmulation.A.data,simmulation.B.data,simmulation.C.data,simmulation.D.data,simmulation.fuel_mass.data,simmulation.ox_mass.data,simmulation.u.Data,simmulation.cg.data(:,1),simmulation.cg.data(:,2),simmulation.cg.data(:,3),simmulation.velocity.data(:,1),simmulation.velocity.data(:,2),simmulation.velocity.data(:,3)])

writetable(csvArray, 'godotData.csv');
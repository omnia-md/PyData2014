import mdtraj as md

traj = md.load("frame0.h5")  # From mdtraj/MDTraj/testing/reference/

indices, phi = md.compute_phi(traj)

plot(phi)
title("Phi Backbone Angle")
xlabel("Timestep")
ylabel("Phi [degrees]")
savefig("./phi.png", bbox_inches=None)

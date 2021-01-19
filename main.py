from LHC import LHC
import matplotlib.pyplot as plt
import numpy as np

machine = LHC(machine_configuration='6.5_TeV_collision_tunes', n_segments=29, D_x=10)

macroparticlenumber_track = 50000
intensity = 1e11
epsn_x  = 2.5e-6
epsn_y  = 3.5e-6
sigma_z = 0.05

bunch   = machine.generate_6D_Gaussian_bunch_matched( macroparticlenumber_track, intensity, epsn_x, epsn_y, sigma_z=sigma_z)

n_turns = 512
beam_x = []
beam_y = []

for i_turn in range(n_turns):
    print('Turn %d/%d'%(i_turn, n_turns))
    machine.track(bunch)

    beam_x.append(bunch.mean_x())
    beam_y.append(bunch.mean_y())


fig, axes = plt.subplots(2, sharex=True)
axes[0].plot(np.array(beam_x))
axes[1].plot(np.array(beam_y))

plt.show()
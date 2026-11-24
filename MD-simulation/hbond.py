import itertools
import mdtraj as md
import numpy as np
import matplotlib.pyplot as plt

#legend_properties = {'weight':'bold',}

traj = md.load('mol20_last_5nstrj.pdb')
#print(t)i

# :func:`md.baker_hubbard` idenfies hydrogen bonds baced on cutoffs
# for the Donor-H...Acceptor distance and angle. The criterion employed
# is :math:`\theta > 120` and :math:`r_\text{H...Acceptor} < 2.5 A` in
# at least 10% of the trajectory. The return value is a list of the 
# indices of the atoms (donor, h, acceptor) that satisfy this criteria.

hbonds = md.baker_hubbard(traj,freq=0.0)
#hbonds = md.wernet_nilsson(traj)
label = lambda hbond : '%s -- %s' % (traj.topology.atom(hbond[0]), traj.topology.atom(hbond[2]))
for hbond in hbonds:
    print(label(hbond))

## Let's compute the actual distances between the donors and acceptors

da_distances = md.compute_distances(traj, hbonds[:, [0,2]], periodic=False)


##print(len(da_distances))
##exit()
#for i in [189]:
#    print(label(hbonds[i]))

#exit()

## Plot a histogram for a few of them
plt.figure()
color = itertools.cycle(['red', 'blue'])
for i in [178 , 189]:
#    #result = [x/1.0 for x in da_distances[:, i]]
#    #print(da_distances[:, i])
#    #exit()
    array = label(hbonds[i]).split("--")
    plt.hist(da_distances[:, i], color=next(color), label=array, alpha=0.5, bins=10,edgecolor='black')

plt.ylabel('Frequency',fontsize=20);
plt.xlabel('Donor-acceptor distance [nm]',fontsize=20)
#plt.ylim(0,35)
#plt.xlim(0.2,0.7)
plt.xticks(fontsize=16)
plt.yticks(fontsize=18)
plt.legend(loc='best',fontsize="large",prop=dict(weight='bold'))
plt.tight_layout()

plt.show()

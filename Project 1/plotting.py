import matplotlib.pyplot as plt 
import numpy as np


def read_file(part):
	filename = "results/part_g/N_10_" + part + ".dat"
	with open(filename, "r") as infile:
		anticipated_max = float(infile.readline())
		tags = infile.readline().strip().split()
		data = {}
		M = len(tags)
		for tag in tags:
			data[tag] = []
		for line in infile:
			line = line.strip().split()
			for	i in range(M):
				data[tags[i]].append(float(line[i]))
		for tag in tags:
			data[tag] = np.array(data[tag])
	return data, anticipated_max

OB, anticipated_max = read_file("OB")
T, anticipated_max = read_file("T")

N = len(OB["x"])
r = np.linspace(-anticipated_max, anticipated_max, N)

cycles_times_particles = 10E8
fontsize = 15
for tag in OB:
	plt.figure(figsize=(8 ,5))
	plt.plot(r, OB[tag]/cycles_times_particles, color="MediumVioletRed", label=r"$\Psi_{OB}$")
	plt.plot(r, T[tag]/cycles_times_particles, color="DarkSlateGray", label=r"$\Psi_{T}$")
	plt.xlabel(r"$%s$" %tag, fontsize=20)
	plt.ylabel(r"$\rho(%s)$" %tag, fontsize=20)
	plt.xticks(fontsize=fontsize)
	plt.yticks(fontsize=fontsize)
	plt.xlim(-anticipated_max, anticipated_max)
	plt.legend(fontsize=20)
	plt.tight_layout()
	plt.savefig("results/part_g/plots/part_g_%s.pdf" %tag)
	plt.show()




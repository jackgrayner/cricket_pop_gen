import pandas as pd
import phlash
import os
import matplotlib.pyplot as plt
import numpy as np
import random

# === Define input VCF and scaffolds
vcf_path = "repeat_filtered.vcf.gz"
scaffolds = ["scaffold_5", "scaffold_9","scaffold_14"]

#get list of all samples 
fam = pd.read_csv("allpops.fam", delim_whitespace=True, header=None, names=["population", "sample"])
popmap = dict(zip(fam["sample"], fam["population"]))
sample_list = fam["sample"].tolist()

#select samples....
samples = [sample for sample, pop in popmap.items() if "Oahu_BYU" in pop]
#select 10
samples=random.sample(samples, 10)

#load each scaffold 
template = ("repeat_filtered_{chrom}.vcf.gz")
chroms = []
for chrom in scaffolds:
    path = os.path.join("./", template.format(chrom=f"{chrom}"))
    chroms.append(
        phlash.contig(path, samples=samples, region=f"{chrom}:10000000-100000000")
    )

results = phlash.fit(data=chroms)

times = np.array([dm.eta.t[1:] for dm in results])
# choose a grid of points at which to evaluate the size history functions
T = np.geomspace(times.min(), times.max(), 1000)
Nes = np.array([dm.eta(T, Ne=True) for dm in results])

plt.figure()
plt.plot(T, np.median(Nes, axis=0))
plt.xscale('log')
plt.yscale('log')
plt.savefig("Ne_Oahu_BYU_chr5_9_14.png")
plt.close()

#try rescaling by Drosophila mutation rates
times = np.array([dm.rescale(2.8e-09).eta.t[1:] for dm in results])
# choose a grid of points at which to evaluate the size history functions
T = np.geomspace(times.min(), times.max(), 1000)
Nes = np.array([dm.eta(T, Ne=True) for dm in results])

plt.figure()
plt.plot(T, np.median(Nes, axis=0))
plt.xscale('log')
plt.yscale('log')
plt.savefig("Ne_Oahu_BYU_chr5_9_14_gens.png")
plt.close()


df = pd.DataFrame(Nes)
df.to_csv("Nes_Oahu_BYU_Chr5_9_14_test2.csv")

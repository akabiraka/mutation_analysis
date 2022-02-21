from cProfile import label
import sys
sys.path.append("../mutation_analysis")

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

df = pd.read_csv("data/merged/merged_ThermoMutDB.csv")

pdb_ids = df["pdb_id"].unique().tolist()
pdb_id = pdb_ids[0]
print("pdb_id: {}".format(pdb_id))
df_of_a_pdb = df[df["pdb_id"]==pdb_id].sort_values(by=["mutation_event", "ddg"])
mutation_events = df_of_a_pdb["mutation_event"].unique().tolist()
print("mutations: {}".format(df_of_a_pdb.shape))
print("n_mutation_events: {}".format(len(mutation_events)))
df_of_a_pdb.to_csv("data/mutations_of_a_pdb.csv", index=False)

X, mutation_sites, means, stds, mins, maxs = [], [], [], [], [], []
for i, mutation_event in enumerate(mutation_events):
    ddgs_df = df_of_a_pdb[df_of_a_pdb["mutation_event"]==mutation_event][["ddg"]]
    if ddgs_df.shape[0]>1:
        X.append(mutation_event)
        mutation_sites.append(int(mutation_event[1:-1]))
        means.append(np.mean(ddgs_df["ddg"]))
        stds.append(np.std(ddgs_df["ddg"]))
        mins.append(np.min(ddgs_df["ddg"]))
        maxs.append(np.max(ddgs_df["ddg"]))
    # if i==5: break

yerr_mins, yerr_maxs = [], []
for i, mean in enumerate(means):
    yerr_mins.append(means[i]-mins[i])
    yerr_maxs.append(maxs[i]-means[i])

# X = [x for _, x in sorted(zip(mutation_sites, X))]
# means = [x for _, x in sorted(zip(mutation_sites, means))]
# stds = [x for _, x in sorted(zip(mutation_sites, stds))]
# mins = [x for _, x in sorted(zip(mutation_sites, mins))]
# maxs = [x for _, x in sorted(zip(mutation_sites, maxs))]
# yerr_mins = [x for _, x in sorted(zip(mutation_sites, yerr_mins))]
# yerr_maxs = [x for _, x in sorted(zip(mutation_sites, yerr_maxs))]

plt.errorbar(X, means, yerr=stds, fmt=".", color="green", ecolor="lightgreen", lw=10, label="Distribution", alpha=0.5)
plt.errorbar(X, means, yerr=[yerr_mins, yerr_maxs], fmt='.', color="green", ecolor='salmon', lw=1, capsize=2, label="Min-max boundary")

for i, mutation_event in enumerate(mutation_events):
    ddgs_df = df_of_a_pdb[df_of_a_pdb["mutation_event"]==mutation_event][["ddg"]]
    if ddgs_df.shape[0]>1:
        for row in ddgs_df.itertuples():
            plt.scatter(mutation_event, row.ddg, color="red", marker=".", alpha=0.5)
    # if i==5: break

plt.xticks(rotation=45)
plt.legend(loc="best")
plt.xlabel("Mutation")
plt.ylabel("Distribution of ddG")
plt.show()
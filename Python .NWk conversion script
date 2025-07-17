import tskit
import os
import random

# ---- CONFIG ----
filepath = os.path.expanduser("~/Documents/SMP3.tree")
output_nwk_path = os.path.expanduser("~/Documents/SMP3.nwk")

inds_per_pop = 150
valid_pops = [1]  # Corrected for SLiM's population IDs starting at 0

# ---- LOAD TREE SEQUENCE ----
ts = tskit.load(filepath)

# ---- SAMPLE INDIVIDUALS ----
sampled_nodes = []
node_labels = {}

for pop_id in valid_pops:
    pop_inds = [ind for ind in ts.individuals() if ind.population == pop_id and len(ind.nodes) > 0]
    print(f"Population {pop_id} has {len(pop_inds)} individuals with nodes.")

    if len(pop_inds) < inds_per_pop:
        print(f"Skipping population {pop_id} (not enough individuals)")
        continue

    sampled_inds = random.sample(pop_inds, inds_per_pop)

    for i, ind in enumerate(sampled_inds):
        node_id = ind.nodes[0]  # Use one node per individual
        sampled_nodes.append(node_id)
        label = f"pop{pop_id}_ind{i}"
        node_labels[node_id] = label

print(f"Sampled {len(sampled_nodes)} nodes in total.")

# ---- SIMPLIFY AND LABEL ----
simplified_ts = ts.simplify(samples=sampled_nodes, keep_unary=True)
tree = simplified_ts.first()

simplified_labels = {}
for n in tree.nodes():
    pop = simplified_ts.node(n).population
    if n in simplified_ts.samples():
        simplified_labels[n] = node_labels.get(n, f"pop{pop}_sample{n}")
    else:
        simplified_labels[n] = f"pop{pop}_node{n}"

# ---- ROOT DISTANCE CALC ----
def path_length_to_root(tree, node):
    length = 0
    while node != tskit.NULL:
        parent = tree.parent(node)
        if parent == tskit.NULL:
            break
        length += tree.branch_length(node)
        node = parent
    return length

print("\nRoot distances for sampled nodes:")
for node in tree.leaves():
    label = simplified_labels.get(node, f"node{node}")
    dist = path_length_to_root(tree, node)
    print(f"{label} → root path length: {dist:.2f}")


# ---- NEWICK EXPORT ----
roots = list(tree.roots)
if len(roots) == 1:
    newick = tree.as_newick(
        root=roots[0],
        node_labels=simplified_labels,
        include_branch_lengths=True
    )
else:
    subtree_newicks = [
        tree.as_newick(
            root=r,
            node_labels=simplified_labels,
            include_branch_lengths=True
        ).rstrip(";")
        for r in roots
    ]
    newick = "(" + ",".join(subtree_newicks) + ");"

with open(output_nwk_path, "w") as f:
    f.write(newick)
print(f"\nSaved Newick tree to {output_nwk_path}")


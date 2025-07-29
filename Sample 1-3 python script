import tskit
import os
import random

# ---- CONFIG ----
filepath = os.path.expanduser("~/Documents/S3.tree")
output_nwk_path = os.path.expanduser("~/Documents/S3.nwk")

valid_pops = [1,2]  # SLiM population IDs start at 0, so 1 is the second
max_samples_per_generation = 5
min_samples_per_generation = 1

# ---- LOAD TREE SEQUENCE ----
ts = tskit.load(filepath)

# ---- GROUP INDIVIDUALS BY GENERATION ----
sampled_nodes = []
node_labels = {}
node_to_population = {}

gen_to_inds = {}
for ind in ts.individuals():
    if ind.population not in valid_pops:
        continue
    if len(ind.nodes) == 0:
        continue
    gen = ind.time
    gen_to_inds.setdefault(gen, []).append(ind)

sorted_gens = sorted(gen_to_inds.keys())
max_time = max(sorted_gens) if sorted_gens else 0

# ---- SAMPLE INDIVIDUALS WITH RECENCY BIAS ----
for gen in sorted_gens:
    inds = gen_to_inds[gen]
    if len(inds) < min_samples_per_generation:
        continue

    recency_weight = 1.0 - (gen / max_time) if max_time > 0 else 1.0
    num_to_sample = round(recency_weight * max_samples_per_generation)
    num_to_sample = max(min_samples_per_generation, min(num_to_sample, len(inds)))

    sampled = random.sample(inds, num_to_sample)
    for i, ind in enumerate(sampled):
        node_id = ind.nodes[0]
        sampled_nodes.append(node_id)
        label = f"gen{int(gen)}_ind{i}"
        node_labels[node_id] = label
        node_to_population[node_id] = ind.population

print(f"✅ Sampled {len(sampled_nodes)} nodes across {len(sorted_gens)} generations with recency bias.")

# ---- SIMPLIFY TREE ----
simplified_ts = ts.simplify(samples=sampled_nodes, keep_unary=True)
tree = simplified_ts.first()

# ---- BUILD MUTATION TYPE MAP ----
mutation_type_map = {}
for site in simplified_ts.sites():
    mut_type = "m2" if int(site.position) == 5000 else "m1"
    for mut in site.mutations:
        mutation_type_map[mut.id] = mut_type

# ---- MAP NODES TO INHERITED MUTATIONS ----
mutations_by_node = {n: [] for n in tree.nodes()}
for mut in simplified_ts.mutations():
    node = mut.node
    mut_type = mutation_type_map.get(mut.id, "unknown")
    mutations_by_node[node].append(mut_type)

# ---- BUILD TIP LABELS WITH POPULATION AND MUTATIONS ----
simplified_labels = {}
for node in tree.leaves():
    inherited_mutations = set()
    cur_node = node
    while cur_node != tskit.NULL:
        inherited_mutations.update(mutations_by_node.get(cur_node, []))
        cur_node = tree.parent(cur_node)

    label_base = node_labels.get(node, f"node{node}")
    pop = node_to_population.get(node, -1)
    pop_label = f"pop{pop}"

    if inherited_mutations:
        muts = "+".join(sorted(inherited_mutations))
        full_label = f"{label_base}|{pop_label}|{muts}"
    else:
        full_label = f"{label_base}|{pop_label}|no_mut"

    simplified_labels[node] = full_label

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

print("\n📏 Root distances for sampled nodes:")
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

print(f"\n✅ Saved Newick tree with mutation- and population-labeled tips to {output_nwk_path}")

import re, sys, glob, os

def parse_newick(s):
    s = s.strip()
    if s.endswith(';'):
        s = s[:-1]
    pos = [0]
    def parse_clade():
        children = []
        if s[pos[0]] == '(':
            pos[0] += 1
            children.append(parse_clade())
            while s[pos[0]] == ',':
                pos[0] += 1
                children.append(parse_clade())
            assert s[pos[0]] == ')', f"expected ) at {pos[0]}: {s[pos[0]-5:pos[0]+5]}"
            pos[0] += 1
        # read label
        m = re.match(r'[^,():;]*', s[pos[0]:])
        label = m.group(0)
        pos[0] += len(label)
        # strip branch length if present
        name = label.split(':')[0]
        return (name, children)
    tree = parse_clade()
    return tree

def get_leaves(node):
    name, children = node
    if not children:
        return {name}
    leaves = set()
    for c in children:
        leaves |= get_leaves(c)
    return leaves

def get_bipartitions(node, all_leaves):
    name, children = node
    bips = set()
    def collect(n):
        nm, ch = n
        if not ch:
            leafset = frozenset([nm])
        else:
            leafset = frozenset()
            for c in ch:
                leafset |= collect(c)
        # only count non-trivial splits (not root, not single leaf, not all leaves)
        if 1 < len(leafset) < len(all_leaves):
            # canonicalize: use the side not containing an arbitrary fixed leaf
            side = leafset
            other = frozenset(all_leaves) - leafset
            canon = min(side, other, key=lambda x: sorted(x))
            bips.add(canon)
        return leafset
    collect(node)
    return bips

def load_trees_from_file(path):
    with open(path) as f:
        content = f.read()
    tree_strs = [t.strip() + ';' for t in content.split(';') if t.strip()]
    return tree_strs

if __name__ == '__main__':
    user_tree_path = sys.argv[1]
    search_dir = sys.argv[2]

    with open(user_tree_path) as f:
        user_str = f.read()
    user_tree = parse_newick(user_str)
    user_leaves = get_leaves(user_tree)
    user_bips = get_bipartitions(user_tree, user_leaves)
    print(f"User tree: {len(user_leaves)} leaves, {len(user_bips)} non-trivial bipartitions")

    for path in sorted(glob.glob(os.path.join(search_dir, '*.nwk'))):
        tree_strs = load_trees_from_file(path)
        for i, ts in enumerate(tree_strs):
            try:
                t = parse_newick(ts)
                leaves = get_leaves(t)
                if leaves != user_leaves:
                    print(f"{os.path.basename(path)} [tree {i}]: leaf set MISMATCH ({len(leaves)} leaves)")
                    continue
                bips = get_bipartitions(t, leaves)
                rf = len(user_bips ^ bips)
                print(f"{os.path.basename(path)} [tree {i}]: RF={rf}  (shared {len(user_bips & bips)}/{len(user_bips)})")
            except Exception as e:
                print(f"{os.path.basename(path)} [tree {i}]: ERROR {e}")

#
# This code computes the GAP group IDs for the component groups of the groups of type N
# enumerated in Section 5. More precisely, we identify the projective image of the intersection
# with SU(3), and the extension by C2.
#
# This code depends on make_groups.sage.

try:
    assert(groups_computed)
except NameError:
    load("make_groups.sage")

def group_id(G):
     gapG = G.gap()
     assert gapG.IsFinite() # Force GAP to do necessary precomputations
     gapH = gapG.Subgroup([diagonal_matrix([z6, z6, z6, ~z6, ~z6, ~z6])])
     return gapG.FactorGroup(gapH).IdGroup()
 
def group_ids_from_list(grouplist):
     return [(name, group_id(G)) for (name, G) in grouplist]

ids = group_ids_from_list(big_list)
groupiddict = dict(ids)
iddict = {}
for (name, groupid) in ids:
     id1 = tuple(groupid)
     name2 = strip_name(name)
     id2 = tuple(groupiddict[name2])
     if (id1, id2) in iddict:
          iddict[(id1,id2)].append(name)
     else:
          iddict[(id1,id2)] = [name]

print("Group IDs computed")
group_ids_computed = True




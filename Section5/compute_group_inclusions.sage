#
# This SageMath code computes the inclusions between the presented groups given in Section 5,
# as used in the proof of Proposition 5.39.
#
# This code also performs the consistency checks described in Remark 5.43 and Remark 5.44.
#
# This code depends on make_groups.sage and compute_group_ids.sage.

try:
    assert group_ids_computed
except NameError:
    load("compute_group_ids.sage")

short_labels = set([i[0] for i in iddict])
for lab in short_labels:
     G = libgap.SmallGroup(lab)
     subs = G.ConjugacyClassesMaximalSubgroups()
     for i in subs:
          H = i.Representative()
          lab2 = tuple(H.IdGroup())
          assert(lab2 in short_labels)

max_names = ["J(B(3, 4; 4))", "J_s(B(1, 12))", "J(B(3, 4))", "J_s(B(3, 4))",
              "J(B(T, 3))", "J_s(B(T, 3))", "J(B(O, 1))", "J(B(O, 2))",
              "J(D(4, 4))", "J(D(6, 6))", "J(E(216))", "J(E(168))"]

found_groups = []
# Iterate over maximal groups.
for name in max_names:
    print("Considering subgroups of {}...".format(name))
    G = groupdict[name]
    gapG = G.gap()
    assert gapG.IsFinite()
    K = gapG.Subgroup([diagonal_matrix([z6, z6, z6, ~z6, ~z6, ~z6])])
    # Pick out subgroups containing K, and in each case identify the isomorphism class of the quotient.
    subs = gapG.ConjugacyClassesSubgroups()
    l = []
    for H0 in subs:
        H = H0.Representative()
        if H.Intersection2(K) != K:
            l.append((H, None))
            continue
        lab1 = H.FactorGroup(K).IdGroup()
        l.append((H, lab1))
    # Identify the index-2 subgroup of G.
    G0 = groupdict[strip_name(name)]
    lab0 = groupiddict[strip_name(name)]
    l1 = [(H, lab1) for (H, lab1) in l if lab1 == lab0 and all(g.matrix() in H for g in G0.gens())]
    assert len(l1) == 1
    G1 = l1[0][0]
    l2 = []
    for (H, lab1) in l:
        if lab1 == None:            
            continue
        lab2 = H.Intersection2(G1).FactorGroup(K).IdGroup()
        # Identify the isomorphism class of H/K and of its index-at-most-2 subgroup.
        l2.append((H, lab1, lab2))
        # Look up all groups matching the two isomorphism classes.
        names = iddict[(tuple(lab1), tuple(lab2))]
        if len(names) > 1:
            lab3 = H.IdGroup()
        for name2 in names:
            # Make sure we didn't find a different maximal subgroup!
            assert not (name2 in max_names and name2 != name)
            # Skip this group if we found it already.
            if name2 in found_groups:
                continue
            # If we have a group that seems to match...
            if len(names) == 1 or groupdict[name2].gap().IdGroup() == lab3:
                # Check whether we have a genuine inclusion.
                if all(g.matrix() in groupdict[name] for g in groupdict[name2].gens()):
                    # If so, mark the subgroup as having been found.
                    found_groups.append(name2)
                else:
                    print("Indirect inclusion of {} into {}".format(name2, name))
      
print("Groups not covered by direct inclusions:")              
print([name for name in groupdict if not name in found_groups])


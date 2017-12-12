Code for 'Robust Vertex Enumeration for Convex Hulls in High Dimensions'
Pranjal Awasthi, Bahman Kalantari, Yikai Zhang
To run experiments, include all folders with subfolders in matlab.
The functions are all included in the '\code' folder.
The AVTA is applied and evaluated in following interesting problems:

Computing vertices: '\computing vertices'
'quick_hull_AVTA.m': AVTA vs Quick hull on computing all vertices problem.
'AVTA_FAW_Simplex.m': AVTA vs Fast Anchor Word on computing vertiecs of perturbed simplex problem.
'AVTA_FAW_Polytope.m': AVTA vs Fast Anchor Word on computing vertiecs of perturbed general convex hull.

Convex hull membership problem: '\convexhull query'
'membership_different_number_pts.m':   Compare the efficiency of AVTA, frank_wolfe, Traignle algorithm and the simplex method on 
solving convexhull query problem. 

Non negative linear system: '\Overcomplete linear system'
'20 query experiments.m': Compare the efficiency of AVTA+ simplex and the simplex method on non negative linear system feasibility problem. 


Non negative matrix factorization(Swimmer dataset): '\non negative matrix factorization\swimmer'
'\swimmer\run_swimmer.m':  Compare AVTA and Fast anchor word on 'Swimmer image ' dataset.


Topic modeling: '\topic modeling'
'\nips\real data\nips_coherence.m': Compute coherence of topics learned by algorithms in 'NIPS' dataset.
'\nips\semi_syhtetic\diff_documents.m': Compute l1 reconstruct error of topics learned by algorithms on NIPS.
'\nips\semi_syhtetic\scaleofnoise.m': Compute the range of l1-error of topics learned by algorithms on NIPS when data has uniform iid noise.
'\nytime\nytime_coherence.m' : Compute coherence of topics learned by algorithms in 'NYtime' dataset.

For bugs and questions please email: yz422@cs.rutgers.edu. Thank you!

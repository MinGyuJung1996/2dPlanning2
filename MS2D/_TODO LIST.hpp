#pragma once
/*
TODO LIST

1. bi1.png
	wrong edges in offCso, offCon
	

2. (solved?)
	wrong edge between connectSlice
	=> guess it was a bug in diff module

3. (SOLVED) : my mistake of thinking (confused a full arc with domain-cut-arcs)
	sn : 0 && i : 658
	cycle is wrong.

4. (done) Support functions
	support(arc)
	support(arcs) := convex hull
	support(model1, model2)
	offset_convexhull

5. max touch circle for non-g1 points
	+) max touch between r > 100.0;

6. (solved) bi2.png
	wrong offset curves
	== was used instead of =

7. (sol?) bi3.png
	mink breaks at 0, 180 degree => guess it is something with the case clause in arc-conv
	temporarilu resolved by giving perturbation to scene(rotation 0.05 deg more)

8. ghost arcs in convex hull

9. voronoi btw cvx Hull
	refine? smaller r?

10. (!) comTan missing.

11. graph construction with new geometry

12. make new scenes

13. (DONE) add comTan to graph(better)
	=> Error in finding line intersection
	=> rotate does not change "Point" itself, but normalize did.

13. dijkstra 
	
*/
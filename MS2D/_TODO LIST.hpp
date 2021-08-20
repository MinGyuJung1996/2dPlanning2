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

14. dijkstra 

15. better lib
	Mink : AS line support (r = inf?)
	Vor: cusp support( need to change transition map into line segments, to support stuff like interval [1.0, 1.0])
	Polymorphism (point/vec) / (arc,line)

16. bi5 ( ve with same arc Idx)
	how?

17. add func: see edges in motion planning
18. fix robot animation
19. how to make better start/end connection?
20. optimize SearchGraph2 (listS?)
21. (DONE) search graph2 : when robot can rotate 360 -> end/first angle connection is needed
22. change elist to take property

23. fig12(edge strange)
24. make more scenes + image
	0. (vor only path) + opt
	1. trail robot
	2. object's with font
	3. read file(model editor)
	4. pick loop(model editor) + delete
	5. add boundary?

25. for now
	path opt more
	more features on scene editor (align global, align diag, zoom-camera, better picking curvatrue)

26. for now
	more straight line -> image refine
		=> display background
			=> zoom/trans cam
	performance experiments
	modelEditor:draw -cspae => change mode
	RSV
	(path-opt) -> AS-path -> rot with arc -> sweep path

27. is arcPath not seeming like it is g1... because some biarcs have eps?
	good way to make polygon biarc

28. completely rework...

29. arcPath rotation at start/end
*/
var documenterSearchIndex = {"docs":
[{"location":"","page":"Home","title":"Home","text":"CurrentModule = GeometricSolutions","category":"page"},{"location":"#GeometricSolutions","page":"Home","title":"GeometricSolutions","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Documentation for GeometricSolutions.","category":"page"},{"location":"","page":"Home","title":"Home","text":"","category":"page"},{"location":"","page":"Home","title":"Home","text":"Modules = [GeometricSolutions]","category":"page"},{"location":"#GeometricSolutions.EnsembleSolution","page":"Home","title":"GeometricSolutions.EnsembleSolution","text":"EnsembleSolution: Collection of all solutions of an EnsembleProblem.\n\n\n\n\n\n","category":"type"},{"location":"#GeometricSolutions.GeometricSolution","page":"Home","title":"GeometricSolutions.GeometricSolution","text":"GeometricSolution: Solution of a geometric differential equation\n\nContains all fields necessary to store the solution of a GeometricProblem.\n\nFields\n\nt:  time steps\ns:  NamedTuple of DataSeries for each solution component\nstep: store every step'th time step (default: 1)\nnstore: number of time steps to store\noffset: offset of current time step\n\nConstructors\n\nGeometricSolution(problem, step = 1)\n\nThe usual way to initialise a Solution is by passing a GeometricProblem, which can for example be an ODEProblem or PODEProblem. The optional parameter step determines the intervals for storing the solution, i.e., if step > 1 only every step'th solution is actually stored.\n\n\n\n\n\n","category":"type"}]
}

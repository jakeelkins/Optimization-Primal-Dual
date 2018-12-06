# Optimization-Primal-Dual
primal-dual algorithm for linear optimization problems


this is my code for the Primal-Dual algorithm for linear optimization problems, a common approach to solve these problems
upon being given an initial "dual feasible" solution--a solution that solves the constraints of program's dual problem,
but not the primal (is "infeasible" for the primal/original problem). 

this algorithm takes matrices A, b, c, and y and returns optimal x and optimal z.

problem: z = c'x s.t. Ax=b, all x >= 0. the initial dual-feasible solution is y.

the algorithm basically creates a tableau from what is called the "associated restricted primal" problem and iteratively
solves that problem to find the final optimal solution (or if the problem is unbouned-->infeasibility).

my .txt files 'a.txt', 'b.txt', 'c.txt', and 'y.txt' are stored as examples for how the program can read more complicated
matrices from a file.

this was turned in for our final project.

Jake Elkins, Math 420 - Linear Optimization Theory, Dr. Min Sun, The University of Alabama, Fall 2018

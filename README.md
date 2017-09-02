# PDC-summer-course

Repository with the project work for a summer course on supercomputing at the
Royal Institute of Technology in Stockholm.

We investigate various ways to speed up a common problem in macroeconomics,
solving the Bellman equation with value function iteration. The problem boils down
to fixed point iteration on a grid, which we do on a Cray supercomputer, taking
race conditions, grid partition (and cache efficiency), etc, into account.

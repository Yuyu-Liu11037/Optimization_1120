import cplex

# ============================================================
# This file uses Cplex Python API to
# establish a Mixed Integer Linear Programming model and then solve it.
# The  problem displayed bellow is as:
#                  min z = cx
#    subject to:      Ax = b
#    and some of x is integer or binary
# ============================================================

# ============================================================
# Input all the data and parameters here

# ============================================================
# Establish the Linear Programming Model
myProblem = cplex.Cplex()

# Add the decision variables and set their lower bound and upper bound
num_decision_var = 105

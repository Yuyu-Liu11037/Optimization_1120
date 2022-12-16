import cplex

# ============================================================
# This file uses Cplex Python API to
# establish an Integer Second Order Cone Programming model and then solve it.
# The  problem displayed bellow is as:
#                  max      z = cx - Omega * sqrt(sigma^2 * x^2)
#           subject to      Ax <= b
#                           and each component of x is binary
#
# The user cuts are
# ============================================================

# ============================================================
# Input all the data and parameters here (随手放了一些数字)
num_decision_var = 26
A = list(range(1, 26))
b = 0.5 * sum(A)
sigma = list(range(1, 26))
c = list(range(1, 26))
Omega = 1

# ============================================================
# Establish the Second Order Cone Programming Model
myProblem = cplex.Cplex()

# Add the decision variables and set their lower bound and upper bound
myProblem.variables.add(names=["x" + str(i) for i in range(num_decision_var - 1)])
myProblem.variables.add(names=["c"])

# Set the type of each variable
for i in range(num_decision_var):
    myProblem.variables.set_types(i, myProblem.variables.type.binary)

# Add constraints
myProblem.linear_constraints.add(
    lin_expr=[cplex.SparsePair(ind=[j for j in range(num_decision_var - 1)], val=A)],
    rhs=[b],
    names=["constraint_1"],
    senses=["L"]
)
myProblem.linear_constraints.add(
    lin_expr=[cplex.SparsePair(ind=26, val=1)],
    rhs=[0],
    names=["constraint_2"],
    senses=["R"]
)

# Add objective function and set its sense
for i in range(num_decision_var):
    myProblem.objective.set_linear([(i, sigma[i])])

myProblem.objective.set_sense(myProblem.objective.sense.maximize)

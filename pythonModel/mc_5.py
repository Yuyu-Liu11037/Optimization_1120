import cplex

# ============================================================
# This file uses Cplex Python API to establish a Mixed Integer Linear
# Programming model and then solve it with User Cuts added in the root node.
# The  problem displayed bellow is as:
#                  min z = cx
#    subject to:      Ax = b
#    and some of x is integer or binary
# ============================================================

# ============================================================
# Input all the data and parameters here
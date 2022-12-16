from docplex.mp.model import Model
from cplex.callbacks import LazyConstraintCallback, UserCutCallback
from docplex.mp.callbacks.cb_mixin import *
import argparse
from docplex.mp.solution import SolveSolution
import sys

class MyCutCallback(ConstraintCallbackMixin, LazyConstraintCallback):
    def __init__(self, env):
        LazyConstraintCallback.__init__(self, env)
        ConstraintCallbackMixin.__init__(self)
        self.eps = 1e-6
        self.nb_cuts = 0

    @property
    def A(self):
        return self.model.A

    @property
    def omega(self):
        return self.model.omega

    # @print_called("--> custom cut callback called: #{0}")
    def __call__(self):
        def f(x, z, I):
            return sum(x[i] * self.A[i] for i in I) - self.omega * z

        # fetch variable solution values at this point.
        s = self.make_complete_solution()

        x_hat = s.get_value_list(self.model.x)
        z_hat = s.get_value(self.model.z)
        t_hat = s.get_value(self.model.obj)

        J = list(range(len(x_hat)))
        J.sort(key=lambda i: x_hat[i], reverse=True)

        h = [f(x_hat, z_hat, J[:n+1]) for n in range(len(x_hat))]

        pi = [b-a for a, b in zip([0]+h, h)]

        if sum(p * x for p, x in zip(pi, x_hat)) > t_hat:
            self.register_constraint(sum(p * x for p, x in zip(pi, self.model.x)) <= t_hat)

        # fetch those constraints which are not satisfied.
        unsats = self.get_cpx_unsatisfied_cts(self.cts, s, self.eps)
        for ct, cut, sense, rhs in unsats:
            # Method add() here is CPLEX's CutCallback.add()
            self.add(cut, sense, rhs)
            self.nb_cuts += 1
            print('-- add new cut[{0}]: [{1!s}]'.format(self.nb_cuts, ct))

def build_model(mdl: Model, n=100, m=5, omega=3):

    A = [i+1 for i in range(n+m)]
    B = [i+1 for i in range(n+m)]
    C = [i+1 for i in range(n+m)]
    b = sum(A) / 2

    x = mdl.binary_var_list(range(1, n+1), name='x_int') + mdl.continuous_var_list(range(n+1, n+m+1), lb=0, ub=1, name='x_rel')
    z = mdl.continuous_var(0, name='z')

    mdl.add(sum(C[i] * x[i] for i in range(n+m)) <= b)
    mdl.add(z * z >= sum(B[i] * B[i] * x[i] * x[i] for i in range(n+m)))
    obj = mdl.sum(A[i] * x[i] for i in range(n+m)) - omega * z
    mdl.maximize(obj)

    mdl.x = x
    mdl.z = z
    mdl.obj = obj
    mdl.A = A
    mdl.omega = omega

    return mdl

def add_callback(mdl: Model):
    return mdl.register_callback(MyCutCallback)

if __name__ == '__main__':
    parser = argparse.ArgumentParser('Model')
    parser.add_argument('-c', '--cut', action='store_true', help='To indicate use cut or not')
    parser.add_argument('-n', type=int, default=25, help='n')
    parser.add_argument('-m', type=int, default=0, help='m')
    parser.add_argument('-e', '--omega', type=int, default=3, help='omega')
    parser.add_argument('-o', '--output', type=str, default='output.txt', help='Output file.')
    parser.add_argument('-w', '--worker', type=int, default=2, help='Number of Threads')
    parser.add_argument('-t', '--timelimit', type=int, default=600, help='Time Limit in seconds.')
    parser.add_argument('-l', '--logout', type=str, default='output.log', help='Output solving log.')

    args = parser.parse_args()

    with Model(name="Model", log_output=True, float_precision=6) as mdl:
        build_model(mdl, n=args.n, m=args.m, omega=args.omega)
        mdl.print_information()
        mdl.parameters.threads = args.worker
        mdl.parameters.timelimit = args.timelimit

        if args.cut:
            print('Use custom cut !')
            cut_cb = add_callback(mdl)
            print('Solving...')

        with open(args.logout, 'w') as f:
            
            stdout_bak = sys.stdout
            sys.stdout = f
            s: SolveSolution = mdl.solve()
            if s:
                x = [v.solution_value for v in mdl.x]
                z = mdl.z.solution_value

                content = [
                    f"Objective value = {mdl.objective_value:5f}",
                    f"x = {x}",
                    f"z = {z}",
                    '===============',
                    "Solver details:",
                    f"status = {s.solve_details.status}",
                    f"gap = {s.solve_details.mip_relative_gap}",
                    f"problem_type = {s.solve_details.problem_type}",
                    f"time = {s.solve_details.time}",
                    f"nb_nodes_processed = {s.solve_details.nb_nodes_processed}",
                    f"columns = {s.solve_details.columns}",
                    f"nb_iterations = {s.solve_details.nb_iterations}"
                    ]

                if args.cut:
                    content.append(f"nb_custom_constraints = {len(cut_cb.cts)}")
                    content.append(f"nb_custom_cuts = {cut_cb.nb_cuts}")

                with open(args.output, 'w') as f:
                    f.write("\n".join(content))

            else:
                print('Model has no solution!')
            
            sys.stdout = stdout_bak

        print(f"\nModel solved, please check '{args.output}' and '{args.logout}'")

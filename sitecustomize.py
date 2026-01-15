import importlib
import importlib.util

spec = importlib.util.find_spec("pulp")
if spec is not None:
    pulp = importlib.import_module("pulp")
    if not hasattr(pulp, "list_solvers") and hasattr(pulp, "listSolvers"):
        pulp.list_solvers = pulp.listSolvers

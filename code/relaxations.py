from gpkit.constraints.relax import ConstantsRelaxed, ConstraintsRelaxedEqually
from gpkit import Model, units
from gpkit.small_scripts import mag

from gpkit.nomials import PosynomialInequality, SignomialInequality
from gpkit.constraints.tight import Tight

import numpy as np

"""
Methods to precondition an SP so that it solves with a relaxed constants algorithm
and postcondition an SP to ensure all relax values are 1
"""

def relaxed_constants(model, include_only=None, exclude=None):
    """
    Method to precondition an SP so it solves with a relaxed constants algorithm

    ARGUMENTS
    ---------
    model: the model to solve with relaxed constants

    RETURNS
    -------
    feas: the input model but with relaxed constants and a new objective
    """

    if model.substitutions:
        constsrelaxed = ConstantsRelaxed(model, include_only, exclude)
        feas = Model(constsrelaxed.relaxvars.prod()**20 * model.cost,
                     constsrelaxed)
        # NOTE: It hasn't yet been seen but might be possible that
        #       the model.cost component above could cause infeasibility
    else:
        feas = Model(model.cost, model)

    return feas

def relaxed_constraints(model):
    """
    Method to precondition an SP so it solves with a relaxed constants algorithm

    ARGUMENTS
    ---------
    model: the model to solve with relaxed constants

    RETURNS
    -------
    feas: the input model but with relaxed constants and a new objective
    """

    constsrelaxed = ConstraintsRelaxedEqually(model)
    return Model(constsrelaxed.relaxvar**20 * model.cost, constsrelaxed)

def aug_dict(key, value, dict):
    if key in dict.keys():
        dict[key] += [value]
    else:
        dict[key] = [value]
    return dict

def post_process(sol):
    """
    Model to print relevant info for a solved model with relaxed constants
    
    ARGUMENTS
    --------
    sol: the solution to the solved model
    """
    print "Checking for relaxed constants..."
    for i in range(len(sol.program.gps)):
        varkeys = [k for k in sol.program.gps[i].varlocs if "Relax" in k.models and sol.program.gps[i].result(k) >= 1.00001]
        if varkeys:
            print "GP iteration %s has relaxed constants" % i
            print sol.program.gps[i].result.table(varkeys)
            if i == len(sol.program.gps) - 1:
                print  "WARNING: The final GP iteration had relaxation values greater than 1"

def compute_constr_tightness(m,sol):
    """
    :param m: model
    :param sol: solution
    :return: dict of inequality constraints and tightness
    Note: only works if all tight constraints have names
    """
    variables = sol['variables']
    tightnessDict = {}
    count = 0
    for constraint in m.flat(constraintsets=True):
        if isinstance(constraint, Tight):
            for i in constraint.flat(constraintsets = False):
                if isinstance(i, PosynomialInequality):
                    leftsubbed = i.left.sub(variables).value
                    rightsubbed = i.right.sub(variables).value
                    rel_diff = abs(1 - leftsubbed/rightsubbed)
                    tightnessDict = aug_dict((constraint.name, count), [rel_diff, i], tightnessDict)
                elif isinstance(i, SignomialInequality):
                    siglt0, = i.unsubbed
                    posy, negy = siglt0.posy_negy()
                    posy = posy.sub(variables).value
                    negy = negy.sub(variables).value
                    rel_diff = abs(1 - posy/negy)
                    tightnessDict = aug_dict((constraint.name, count), [rel_diff, i], tightnessDict)
            count += 1
    return tightnessDict

def group_constr_tightness(tightnessDict):
    strkeys = list(set([key[0] for key in tightnessDict.keys()]))
    groupedDict = {}
    for key, value in sorted(tightnessDict.iteritems()):
        if key[0] in groupedDict.keys():
            groupedDict[key[0]] += value
        else:
            groupedDict[key[0]] = value
    return groupedDict

def map_relaxations(groupedDict, nt, nx, decimals = 3):
    strkeys = groupedDict.keys()
    relaxDict = {i:None for i in strkeys}
    for i in strkeys:
        dim = len(groupedDict[i])
        nd = groupedDict[i]
        # Mapping relaxations in time and space
        dim1, dim2 = [0,0]
        if dim == nx*nt:
            [dim1,dim2] = [nt,nx]
        elif dim == nt:
            [dim1,dim2] = [nt, 1]
        elif dim == nx:
            [dim1,dim2] = [1, nx]
        # Special cases where monomials (always tight) replace posys
        # in the first section
        elif dim == (nx-1)*nt and i == 'massCons':
            [dim1,dim2] = [nt, nx]
            a = np.array(nd)[:,0].reshape((nt, nx-1))
            nd = np.concatenate((np.zeros((nt,1))*units(''), a), axis=1)
            nd = nd.reshape((nx*nt,1))
        else:
            print 'Warning: tight constraint that does not' \
                  ' obey the specified relaxations detected.'
        nd = [mag(j) for j in np.array(nd)[:,0]]
        nd = np.array(nd).reshape((dim2, dim1))
        relaxDict[i] = np.round(nd, decimals=decimals)
    return relaxDict








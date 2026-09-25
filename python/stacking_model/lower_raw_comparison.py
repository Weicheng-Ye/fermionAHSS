"""Explicit legal beta versus natural raw-beta comparison in input degree3.

The signed-closed-A lower bridge extends the mismatch to a closed binary
operation on arbitrary B. Its A-only source uses the existing shared-sign
pair contraction. The relative degree5 complex is acyclic modulo2, so no
source coefficients or cochain equations are fitted on X.
"""
from functools import lru_cache
from itertools import combinations
import off_shell_beta as off
import closed_a_upper as sharp
import nonzero_offshell_comparison as bridge
import stacking_lower as lower
import v3_pair_shared as model
from g6_repair import f2_solve
p,hp,rc=off.p,off.hp,model.rc


def extended_beta(A,B,Ap,Bp,s,omega):
    """Natural beta for fhat on signed-closed A and arbitrary B."""
    alpha=off.alpha(A,Ap,s)
    As,Bs=A+Ap,p.binary(B+Bp+alpha)
    b=p.binary(p.differential(B)+off.primary(A,s,omega))
    bp=p.binary(p.differential(Bp)+off.primary(Ap,s,omega))
    return p.binary(sharp.beta_sharp(A,B,Ap,Bp,s,omega)
                    +bridge.lower_gauge(A,B,s,omega)
                    +bridge.lower_gauge(Ap,Bp,s,omega)
                    +bridge.lower_gauge(As,Bs,s,omega)
                    +p.cup(b,bp,b.degree))


def residual(A,B,Ap,Bp,s,omega):
    """Closed degree-(n+2) binary cochain for signed-closed A,A'."""
    return p.binary(extended_beta(A,B,Ap,Bp,s,omega)
                    +off.raw_beta(A,B,Ap,Bp,s,omega,lower.legal_beta))


def source(A,Ap,s,omega):
    return residual(A,p.zero(A.degree+1),Ap,p.zero(A.degree+1),s,omega)


@lru_cache(None)
def relative_acyclicity():
    low,mid,high=(rc.total_basis(q)for q in(4,5,6))
    if low or len(mid)!=2 or len(high)!=12:
        raise ArithmeticError('the prescribed relative degree5 complex changed')
    columns=tuple(sum((rc.total_boundary(key).get(cell,0)%2)<<j
                      for j,cell in enumerate(mid))for key in high)
    rows=tuple(sum(((column>>j)&1)<<i for i,column in enumerate(columns))
               for j in range(len(mid)))
    witnesses=tuple(f2_solve(rows,1<<j,len(high))for j in range(len(mid)))
    if any(witness is None for witness in witnesses):
        raise ArithmeticError('the relative binary degree5 source is not acyclic')
    return dict(lower=low,middle=mid,upper=high,boundary_columns=columns,
                cycle_fillings=witnesses)


@lru_cache(None)
def source_value(pair):
    A,Ap,s=model.from_diag(pair[0]);vertices=tuple(range(6))
    if (not any(A(f)for f in combinations(vertices,4))
            or not any(Ap(f)for f in combinations(vertices,4))):
        return 0
    return source(A,Ap,s,model.from_omega(pair[1]))(vertices)


def evaluate(chain):
    return sum(c*source_value(pair)for pair,c in chain.items())%2


@lru_cache(None)
def small_source_values():
    data=relative_acyclicity()
    values=tuple(evaluate(rc.Gtot(cell))for cell in data['middle'])
    if any(values):
        raise ArithmeticError('the closed normalized comparison violated its acyclic source certificate')
    return values


def A_primitive(A,Ap,s,omega):
    if (A.degree,Ap.degree,s.degree,omega.degree)!=(3,3,1,2):
        raise ValueError('this fixed comparison constructor expects degrees(3,3,1,2)')
    small_source_values()
    def value(vertices):
        # The fixed contraction preserves both axes, and R_A is literally
        # zero there. Avoid constructing a large zero-valued chain.
        if (not any(A(face)for face in combinations(vertices,4))
                or not any(Ap(face)for face in combinations(vertices,4))):
            return 0
        pair=(model.to_diag(A,Ap,s,vertices),model.to_omega(omega,vertices))
        return evaluate(rc.Htot(pair))
    return p.Cochain(4,value)


def primitive(A,B,Ap,Bp,s,omega):
    """lambda3 with delta lambda3=beta_legal+beta_raw on legal AB pairs."""
    I=hp.interval
    along=residual(I(A),I(B,True),I(Ap),I(Bp,True),I(s),I(omega))
    return p.binary(A_primitive(A,Ap,s,omega)+hp.prism(along))


def local_extension(A,B,Ap,Bp,s,omega):
    """Natural degree-four extension, testing legality on each simplex."""
    if (A.degree,B.degree,Ap.degree,Bp.degree)!=(3,4,3,4):
        raise ValueError('the local comparison expects two degree-(3,4) pairs')
    curvatures=off.curvature(A,B,s,omega)+off.curvature(Ap,Bp,s,omega)
    gauge=primitive(A,B,Ap,Bp,s,omega)
    def value(vertices):
        if any(c(face)for c in curvatures
               for face in combinations(vertices,c.degree+1)):
            return 0
        return gauge(vertices)
    return p.Cochain(4,value)


def natural_beta(A,B,Ap,Bp,s,omega):
    """Natural beta6 equal to the calibrated beta on every legal pair."""
    return p.binary(off.raw_beta(A,B,Ap,Bp,s,omega,lower.legal_beta)
                    +p.differential(local_extension(A,B,Ap,Bp,s,omega)))

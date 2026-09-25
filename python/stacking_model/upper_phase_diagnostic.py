"""Bounded diagnostics only; no asserted new upper stacking formula."""
from functools import lru_cache
import h_tau_primitive as hp
p,low=hp.p,hp.low


def production_phase(n,A,B,C,s,omega):
    def lift(c):
        def value(face):
            base=tuple(v[0] for v in face)
            return 0 if any(a==b for a,b in zip(base,base[1:])) else c(base)
        return p.Cochain(c.degree,value)
    args=tuple(map(lift,(A,B,C,s,omega)))
    if n<3:
        operation=low.build_phase(n,*args)
    else:
        import high_phase
        operation=high_phase.build_phase(*args)
    return p.Cochain(n+4,lambda f:operation(tuple((v,) for v in f)))


from closed_a_upper import fsharp,beta_sharp


if __name__=='__main__':
    import random,time
    from itertools import combinations
    n=1;z=tuple(range(n+6));r=random.Random(7124)
    eps={i:r.randrange(2)for i in z};s=p.Cochain(1,lambda f:(eps[f[0]]+eps[f[1]])%2)
    wt={f:r.randrange(2)for f in combinations(z,2)};w=p.binary(p.differential(p.Cochain(1,lambda f:wt.get(f,0))))
    at={f:r.randrange(3)-1 for f in combinations(z,n)};da=p.differential(p.Cochain(n-1,lambda f:at.get(f,0)))
    A=p.Cochain(n,lambda f:(-1)**eps[f[0]]*da(f))
    base=production_phase(n,A,p.zero(n+1),p.zero(n+2),s,w)
    src=hp.uniform_source(A,s,w)
    target=p.theta(src['p'],src['k0'],s,w)
    print('start',flush=True);t=time.time()
    print('base derivative',p.ds(base,s)(z),'source',target(z),'sum',(p.ds(base,s)(z)+target(z))%1,'time',time.time()-t,flush=True)

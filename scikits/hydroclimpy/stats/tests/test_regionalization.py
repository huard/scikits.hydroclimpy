"""
Tests for regionalization
Testing against the R library would be good but the lmomco library does not include regional homogeneity functions.
Instead reglmr and regstat are tested against output results incuded in the source distribution of the fortran code. 
"""
# David Huard, 2009

from numpy.ma.testutils import *
from scikits.hydroclimpy.stats import reglmr, regstat, lmoments, discordance
import scipy.stats.distributions as dist

class TestRegionalization(TestCase):
    def test_reglmr_simple(self):
        rmom = reglmr(np.array([[2., 1., .5]]).T, weight=[1],  nmom=3)
        assert_equal(rmom, [1.,  .5,  .5])
        
    def test_reglmr(self):
        nsite = 10 
        L1 = dist.gamma.rvs(4, size=nsite)
        L2 = dist.gamma.rvs(2., size=nsite)
        R3 = dist.beta.rvs(3., 5., size=nsite)
        R4 = dist.beta.rvs(4., 5., size=nsite)
        W = dist.poisson.rvs(4, loc=10, size=nsite)
        xmom = np.vstack([L1, L2, R3, R4])
        rmom = reglmr(xmom, W)
        assert_equal(len(xmom), 4)
        assert_equal(rmom[2], np.average(R3, weights=W))
        
        
    def test_reglmr_maxwind(self):
        """Compare results with those in MAXWIND.OUT obtained by 
        processing the MAXWIND.DAT input file."""
        import maxwind
        
        lmom = {}
        L = {}
        # Compute sample L-moments
        for name, d in maxwind.data.items():
            lmom[name] = lmoments(d, nmom=5, a=-0.35, b=0, mode='stats')
            L[name] = len(d)
            assert_almost_equal(lmom[name], maxwind.lmr[name],2)
            
        # Compute regionalized value.
        xmom = np.array(lmom.values()).T
        w = np.array(L.values())
        assert_almost_equal(np.around(reglmr(xmom, w, 5), 4), maxwind.reglmr,6)
        
    def test_regstat_cascades(self):
        """Compare results with those in CASCADES.OUT obtained by
        processing the CASCADES.DAT input file."""
        import cascades as cas      
        N = cas.input['N']
        lmom = cas.input['lmom'].T
        prob = [.01, .02, .05, .1, .2, .5, .9, .95, .99, .999]
        R = regstat(N, lmom, prob, nsim=500, seed=619145091.)
        assert_almost_equal(R['d'], cas.discordance, 2)
        assert_almost_equal(R['rmom'][1:4], cas.rmom, 2)
        assert_almost_equal(R['vobs'], cas.vobs, 2)
        assert_almost_equal(R['vbar'], cas.vbar, 2)
        assert_almost_equal(R['vobs'], cas.vobs, 2)
        assert_almost_equal(R['h'], cas.h, 2)
        assert_almost_equal(R['para'][:3,2], cas.para_gen_normal, 2)
        assert_almost_equal(R['para'][:3,3], cas.para_pearson, 2)
        assert_almost_equal(R['para'][:,5], cas.para_wakeby, 2)
        
    def test_discordance(self):
        import cascades as cas
        D = discordance(cas.input['lmom'])
        assert_almost_equal(D, cas.discordance, 2)
        
if __name__ == '__main__':
    run_module_suite()
    

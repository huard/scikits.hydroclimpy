##------------------------------------------------------------------------------
## Regional Analysis
##------------------------------------------------------------------------------

import numpy as np
from scikits.hydroclimpy.stats import _lmoments

def reglmr(xmom, weight, nmom=4):
    """Regional weighted average of L-ratios.
    
    Parameters
    ----------
    xmom : array of dimension (nxmom,nsite). 
        xmom(i,j) contains the i'th L-statistics for site j (L-1, L-2, T-3, T-4, ...). 
    weight : array of length (nsite)
        The weights to be applied to each site, typically the 
        length of the series at each site. 
    nmom : int
        Number of L-ratios to be found (default=4).
    
    Retuns
    ------
    rmom   : array of length nmom. 
        The regional weighted average L-ratios (1, T, T-3, T-4, ...), where T is L-2/L-1.
    """
    import numpy as np
    rmom = np.empty(nmom)
    xmom = np.atleast_2d(xmom)
    if xmom.shape[1] != len(weight):
        raise ValueError("Bad shape for xmom, got %s and expected %d in the last dimension."%(str(xmom.shape), len(weight)))
    rmom[0] = 1.
    rmom[1] = np.average(xmom[1, :]/xmom[0, :], weights=weight)
    rmom[2:] = np.average(xmom[2:nmom, :], axis=1, weights=weight)
    return rmom


def discordance(xmom):
    """Return the discordance measure at each site.
    
    Parameters
    ----------
    xmom : ndarray (N, 4)
      First 4 sample L-moments for each site, in the order: mean, L-cv, L-skewness,
      L-kurtosis. Note that L-cv is required, not L-2. 
      
    Returns
    -------
    D : array (N,)
      Discordance measure. Large values indicate potential errors at that site. 
    """
    # We are not using T5. Do we still  ask for it ?
    nsites = xmom.shape[0]
    
    if nsite <= 3:
        raise ValueError("More than three sites are needed, %d given." % nsites)
        
    if xmom.shape[1] != 4:
        raise ValueError("The first four L-moments are required, %d given." % xmom.shape[1])
    
    # Deviation from the mean
    anomaly = np.ma.anom(xmom[:,1:-1], axis=0)
    
    # Sum of squared differences
    smat = np.dot(anomaly.T, anomaly)
    
    # Matrix inverse
    ismat = np.linalg.inv(smat)
    
    # Discordance
    D = np.zeros(nsites)
    for i in range(nsites):
        D[i] = np.dot(np.dot(anomaly[i], ismat), anomaly[i])
    
    return  D * nsites / 3.


def regstat(length,xmom,prob,nsim=1,a=0,b=0,seed=123,kprint=1,kout=1):
    """calculates three statistics useful in regional frequency analysis

 Discordancy measure, d(i), for individual sites in a region.
    Large values might be used as a flag to indicate potential errors
    in the data at the site.  "Large" might be 3 for regions with 15
    or more sites, but less (exact values in array dc1) for smaller
    regions.

 Heterogeneity measures, h(j), for a region based upon either:-
    j=1: the weighted s.d. of the l-cvs or
    j=2: the average distance from the site to the regional average
         on a graph of l-cv vs. l-skewness
    j=3: the average distance from the site to the regional average
         on a graph of l-skewness vs. l-kurtosis
    In practice h(1) is probably sufficient.  A value greater than
    (say) 1.0 suggests that further subdivision of the region should
    be considered as it might improve quantile estimates.

 Goodness-of-fit measures, z(k), for 5 candidate distributions:
    k=1: generalized logistic
    k=2: generalized extreme value
    k=3: generalized normal (lognormal)
    k=4: pearson type iii (3-parameter gamma)
    k=5: generalized pareto
    Provided that the region is acceptably close to homogeneous,
    the fit may be judged acceptable at 10% significance level
    if z(k) is less than 1.645 in absolute value.

 For further details see J.R.M. Hosking and J.R. Wallis (1997),
 "Regional frequency analysis: an approach based on l-moments",
 cambridge university press, chapters 3-5.

 parameters of routine:
    length : array of length nsites. record lengths at each site.
    xmom   : Array of dimension (5,nsites), containing the first 5 sample
             L-moments for each site, in the order mean, L-cv, L-skewness,
             L-kurtosis, T-5, i.e xmom(i,j) contains the i'th l-moment for site j.
             NB xmom(2,.) contains l-cv, not the usual l-2!
    a, b   : Parameters of plotting position.
             Note: a and b should be the same as the values used
             to calculate the moments in the xmom array.
    seed   : Seed for random number generator. Should be a whole number
             in the range 2d0 to 2147483647d0.
    nsim   : Number of simulated worlds for heterogeneity and
             goodness-of-fit tests.
             Note: nsim=0 will force return at completion of
                 outlier test.  nsim=1 will suppress calculation of
                 h and z statistics, but parameter and quantile
                 estimates will be found.
    Prob   : Array of length nprob.  probabilities for which
             quantiles are to be calculated.
    kprint : Output flag. should be set to
                 0  to suppress output
                 1  to print output
    kout   : Channel to which output is directed
Outputs:    
    rmom   : Array of length 5. On exit, contains the regional
             weighted average L-moment ratios.
    d      : Array of length nsites. on exit, contains the
                 discordancy measure (d statistic) for each site.
    vobs   : Array of length 3. On exit, contains the regional
             observed values of 3 heterogeneity statistics:
                 (1) weighted s.d. of L-cvs;
                 (2) average of L-cv/L-skew distances;
                 (3) average of L-skew/L-kurtosis distances.
    vbar   : Array of length 3. On exit, contains the mean of the
             simulated values of the 3 heterogeneity statistics.
    vsd    : Array of length 3. On exit, contains the s.d. of the
             simulated values of the 3 heterogeneity statistics.
    h      : Array of length 3. On exit, contains heterogeneity
             measures (h statistics), i.e. h=(vobs-vbar)/vsd.
    z      : Array of length 5. On exit, contains goodness-of-fit
             measures (z statistics) for 5 distributions:
                 (1) gen. logistic, (2) gen. extreme value,
                 (3) gen. normal, (4) pearson type iii,
                 (5) gen. pareto.
     para   : Array of dimension (5,6). On exit, if nsim.ge.1,
              contains parameters of growth curves fitted by the
              above 5 distributions, plus wakeby.
    """
    #(rmom,d,vobs,vbar,vsd,h,z,para) = 
    out = _lmoments.regtst(length,xmom,a,b,seed,nsim,prob,kprint,kout)
    ierr = out[-1]
    if ierr == -1:
        raise ValueError("Invalid parameter nmom in samlmr.")
    elif ierr == -2:
        raise ValueError("Invalid plotting-position parameters.")
    elif ierr == -3:
        warnings.warn("All data input are equal.")
    elif ierr == -10:
        raise ValueError("Number of sites exceed maximum. Recompile library with MAXNS set to a higher value.")
    elif ierr == -20:
        raise ValueError("Number of quantiles exceed maximum. Recompile library with MAXQ set to a higher value or reduce the length of prob.")
    elif ierr == -30:
        raise ValueError("Number of records exceed maximum. Recompile library with MAXREC set to a higher value.")
    elif ierr == -40:
        warnings.warn("Unable to invert sum-of-squares matrix, the discordance statistics will be computed.")
    
    return dict( zip(['rmom', 'd', 'vobs', 'vbar', 'vsd', 'h', 'z', 'para'], out[:-1]) )
    #return (rmom, d,
    #        np.reshape(np.concatenate((vobs,vbar,vsd)),(3,-1)),
    #        np.reshape(np.concatenate((h,z,np.ravel(para))),(5,8)))
    

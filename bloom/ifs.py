"""
IFS Module
===================
Placeholder.

Authors: Eduardo G. Gusmao.

"""

###################################################################################################
# Libraries
###################################################################################################

# Python

# Internal

# External

###################################################################################################
# Ifs Class
###################################################################################################

class Ifs():
  """This class represents TODO.

  *Keyword arguments:*

    - argument1 -- Short description. This argument represents a long description. It can be:
      - Possibility 1: A possibility 1.
      - Possibility 2: A possibility 2.

    - argument2 -- Short description. This argument represents a long description. It can be:
      - Possibility 1: A possibility 1.
      - Possibility 2: A possibility 2.
  """

  def __init__(self):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """
    pass

  def calculate_ifs(self):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """

    # Create a vector of significant contact points according to the points found in the dicionaries of SICA. Assign a value based on the point and their neighborhood - given the standardized matrix from DPMM.

  def fix_matrix(self):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """

    # 1. Remove intervals outside the chrom dict boundaries
    # 2. Put all negative values to 0 and all >1 to 1
    # 3. Make all sica.removed_dict rows/columns = 0












class Prior():
  """This class represents TODO.

  *Keyword arguments:*

    - argument1 -- Short description. This argument represents a long description. It can be:
      - Possibility 1: A possibility 1.
      - Possibility 2: A possibility 2.

    - argument2 -- Short description. This argument represents a long description. It can be:
      - Possibility 1: A possibility 1.
      - Possibility 2: A possibility 2.
  """

  def __init__(self, post=None, *args, **kwargs):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """
    if post is None:
      post = type(self)
    self._post = post

  def sample(self, size=None):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """
    raise NotImplementedError

  def like1(self, x, *args, **kwargs):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """
    raise NotImplementedError

  def likelihood(self, D, *args, **kwargs):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """
    return np.prod(self.like1(D, *args, **kwargs))

  def lnlikelihood(self, D, *args, **kwargs):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """
    return np.log(self.likelihood(D, *args, **kwargs))

  def __call__(self, *args):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """
    raise NotImplementedError

  def _post_params(self, D):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """
    raise NotImplementedError

  def post(self, D):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """
    return self._post(*self._post_params(D))

  def pred(self, x):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """
    raise NotImplementedError



class GaussianMeanKnownVariance(Prior):
  """This class represents TODO.

  *Keyword arguments:*

    - argument1 -- Short description. This argument represents a long description. It can be:
      - Possibility 1: A possibility 1.
      - Possibility 2: A possibility 2.

    - argument2 -- Short description. This argument represents a long description. It can be:
      - Possibility 1: A possibility 1.
      - Possibility 2: A possibility 2.
  """

  def __init__(self, mu_0, sigsqr_0, sigsqr):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """
    self.mu_0 = mu_0
    self.sigsqr_0 = sigsqr_0
    self.sigsqr = sigsqr
    self._norm1 = np.sqrt(2*np.pi*self.sigsqr)
    self._norm2 = np.sqrt(2*np.pi*self.sigsqr_0)
    super(GaussianMeanKnownVariance, self).__init__()

  def sample(self, size=None):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """
    if size is None:
      return np.random.normal(self.mu_0, np.sqrt(self.sigsqr_0))
    else:
      return np.random.normal(self.mu_0, np.sqrt(self.sigsqr_0), size=size)

  def like1(self, x, mu):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """
    return np.exp(-0.5*(x-mu)**2/self.sigsqr) / self._norm1

  def __call__(self, mu):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """
    return np.exp(-0.5*(mu-self.mu_0)**2/self.sigsqr_0) / self._norm2

  def _post_params(self, D):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """
    try:
      n = len(D)
    except TypeError:
      n = 1
    Dbar = np.mean(D)
    sigsqr_n = 1./(n/self.sigsqr + 1./self.sigsqr_0)
    mu_n = sigsqr_n * (self.mu_0/self.sigsqr_0 + n*Dbar/self.sigsqr)
    return mu_n, sigsqr_n, self.sigsqr

  def pred(self, x):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """
    sigsqr = self.sigsqr + self.sigsqr_0
    return np.exp(-0.5*(x-self.mu_0)**2/sigsqr) / np.sqrt(2*np.pi*sigsqr)



class InvGamma(Prior):
  """This class represents TODO.

  *Keyword arguments:*

    - argument1 -- Short description. This argument represents a long description. It can be:
      - Possibility 1: A possibility 1.
      - Possibility 2: A possibility 2.

    - argument2 -- Short description. This argument represents a long description. It can be:
      - Possibility 1: A possibility 1.
      - Possibility 2: A possibility 2.
  """

  def __init__(self, alpha, beta, mu):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """
    self.alpha = alpha
    self.beta = beta
    self.mu = mu
    super(InvGamma, self).__init__()

  def sample(self, size=None):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """
    return 1./np.random.gamma(self.alpha, scale=self.beta, size=size)

  def like1(self, x, var):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """
    return np.exp(-0.5*(x-self.mu)**2/var) / np.sqrt(2*np.pi*var)

  def __call__(self, var):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """
    al, be = self.alpha, self.beta
    return be**(-al)/gamma(al) * var**(-1.-al) * np.exp(-1./(be*var))

  def _post_params(self, D):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """
    try:
      n = len(D)
    except TypeError:
      n = 1
    al_n = self.alpha + n/2.0
    be_n = 1./(1./self.beta + 0.5*np.sum((np.array(D)-self.mu)**2))
    return al_n, be_n, self.mu

  def pred(self, x):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """
    return t_density(2*self.alpha, self.mu, 1./self.beta/self.alpha, x)

  def evidence(self, D):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """
    raise NotImplementedError



class InvGamma2D(Prior):
  """This class represents TODO.

  *Keyword arguments:*

    - argument1 -- Short description. This argument represents a long description. It can be:
      - Possibility 1: A possibility 1.
      - Possibility 2: A possibility 2.

    - argument2 -- Short description. This argument represents a long description. It can be:
      - Possibility 1: A possibility 1.
      - Possibility 2: A possibility 2.
  """

  def __init__(self, alpha, beta, mu):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """
    self.alpha = alpha
    self.beta = beta
    self.mu = np.array(mu)
    assert len(mu) == 2
    super(InvGamma2D, self).__init__()

  def sample(self, size=None):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """
    return 1./np.random.gamma(self.alpha, scale=self.beta, size=size)

  def like1(self, x, var):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """
    assert isinstance(x, np.ndarray)
    assert x.shape[-1] == 2
    return np.exp(-0.5*np.sum((x-self.mu)**2, axis=-1)/var) / (2*np.pi*var)

  def lnlikelihood(self, D, var):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """
    return -0.5*np.sum((D-self.mu)**2)/var - D.shape[0]*np.log(2*np.pi*var)

  def __call__(self, var):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """
    al, be = self.alpha, self.beta
    return be**(-al)/gamma(al) * var**(-1.-al) * np.exp(-1./(be*var))

  def _post_params(self, D):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """
    try:
      n = len(D)
    except TypeError:
      n = 1
    al_n = self.alpha + n
    be_n = 1./(1./self.beta + 0.5*np.sum((np.array(D)-self.mu)**2))
    return al_n, be_n, self.mu

  def pred(self, x):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """
    assert isinstance(x, np.ndarray)
    assert x.shape[-1] == 2
    return multivariate_t_density(2*self.alpha, self.mu, 1./self.beta/self.alpha*np.eye(2), x)

  def evidence(self, D):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """
    raise NotImplementedError



class NormInvChi2(Prior):
  """This class represents TODO.

  *Keyword arguments:*

    - argument1 -- Short description. This argument represents a long description. It can be:
      - Possibility 1: A possibility 1.
      - Possibility 2: A possibility 2.

    - argument2 -- Short description. This argument represents a long description. It can be:
      - Possibility 1: A possibility 1.
      - Possibility 2: A possibility 2.
  """

  def __init__(self, mu_0, kappa_0, sigsqr_0, nu_0):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """
    self.mu_0 = float(mu_0)
    self.kappa_0 = float(kappa_0)
    self.sigsqr_0 = float(sigsqr_0)
    self.nu_0 = float(nu_0)
    super(NormInvChi2, self).__init__()
    model_dtype = np.dtype([('mu', float), ('var', float)])

  def sample(self, size=None):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """
    if size is None:
      var = 1./np.random.chisquare(df=self.nu_0)*self.nu_0*self.sigsqr_0
      ret = np.zeros(1, dtype=self.model_dtype)
      ret['mu'] = np.random.normal(self.mu_0, np.sqrt(var/self.kappa_0))
      ret['var'] = var
      return ret[0]
    else:
      var = 1./np.random.chisquare(df=self.nu_0, size=size)*self.nu_0*self.sigsqr_0
      ret = np.zeros(size, dtype=self.model_dtype)
      ret['mu'] = (np.random.normal(self.mu_0, np.sqrt(1./self.kappa_0), size=size) * np.sqrt(var))
      ret['var'] = var
      return ret

  def like1(self, *args):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """
    if len(args) == 3:
      x, mu, var = args
    elif len(args) == 2:
      x, theta = args
      mu = theta['mu']
      var = theta['var']
    return np.exp(-0.5*(x-mu)**2/var) / np.sqrt(2*np.pi*var)

  def __call__(self, *args):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """
    if len(args) == 2:
      mu, var = args
    elif len(args) == 1:
      mu = args[0]['mu']
      var = args[0]['var']
    return (normal_density(self.mu_0, var/self.kappa_0, mu) * scaled_IX_density(self.nu_0, self.sigsqr_0, var))

  def _post_params(self, D):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """
    try:
      n = len(D)
    except TypeError:
      n = 1
    Dbar = np.mean(D)
    kappa_n = self.kappa_0 + n
    mu_n = (self.kappa_0*self.mu_0 + n*Dbar)/kappa_n
    nu_n = self.nu_0 + n
    sigsqr_n = ((self.nu_0*self.sigsqr_0 + np.sum((D-Dbar)**2) + n*self.kappa_0/(self.kappa_0+n)*(self.mu_0-Dbar)**2)/nu_n)
    return mu_n, kappa_n, sigsqr_n, nu_n

  def pred(self, x):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """
    return t_density(self.nu_0, self.mu_0, (1.+self.kappa_0)*self.sigsqr_0/self.kappa_0, x)

  def evidence(self, D):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """
    mu_n, kappa_n, sigsqr_n, nu_n = self._post_params(D)
    try:
      n = len(D)
    except:
      n = 1
    return (gamma(nu_n/2.0)/gamma(self.nu_0/2.0) * np.sqrt(self.kappa_0/kappa_n) *
           (self.nu_0*self.sigsqr_0)**(self.nu_0/2.0) / (nu_n*sigsqr_n)**(nu_n/2.0) / np.pi**(n/2.0))

  def marginal_var(self, var):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """
    return scaled_IX_density(self.nu_0, self.sigsqr_0, var)

  def marginal_mu(self, mu):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """
    return t_density(self.nu_0, self.mu_0, self.sigsqr_0/self.kappa_0, mu)



class NormInvGamma(Prior):
  """This class represents TODO.

  *Keyword arguments:*

    - argument1 -- Short description. This argument represents a long description. It can be:
      - Possibility 1: A possibility 1.
      - Possibility 2: A possibility 2.

    - argument2 -- Short description. This argument represents a long description. It can be:
      - Possibility 1: A possibility 1.
      - Possibility 2: A possibility 2.
  """

  def __init__(self, m_0, V_0, a_0, b_0):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """
    self.m_0 = float(m_0)
    self.V_0 = float(V_0)
    self.a_0 = float(a_0)
    self.b_0 = float(b_0)
    super(NormInvGamma, self).__init__()
    model_dtype = np.dtype([('mu', float), ('var', float)])

  def sample(self, size=None):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """
    if size is None:
      var = 1./np.random.gamma(self.a_0, scale=1./self.b_0)
      ret = np.zeros(1, dtype=self.model_dtype)
      ret['mu'] = np.random.normal(self.m_0, np.sqrt(self.V_0*var))
      ret['var'] = var
      return ret[0]
    else:
      var = 1./np.random.gamma(self.a_0, scale=1./self.b_0, size=size)
      ret = np.zeros(size, dtype=self.model_dtype)
      ret['mu'] = np.random.normal(self.m_0, np.sqrt(self.V_0), size=size)*np.sqrt(var)
      ret['var'] = var
      return ret

  def like1(self, *args):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """
    if len(args) == 3:
      x, mu, var = args
    elif len(args) == 2:
      x, theta = args
      mu = theta['mu']
      var = theta['var']
    return np.exp(-0.5*(x-mu)**2/var) / np.sqrt(2*np.pi*var)

  def __call__(self, *args):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """
    if len(args) == 1:
      mu = args[0]['mu']
      var = args[0]['var']
    elif len(args) == 2:
      mu, var = args
    normal = np.exp(-0.5*(self.m_0-mu)**2/(var*self.V_0))/np.sqrt(2*np.pi*var*self.V_0)
    ig = self.b_0**self.a_0/gamma(self.a_0)*var**(-(self.a_0+1))*np.exp(-self.b_0/var)
    return normal*ig

  def _post_params(self, D):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """
    try:
      n = len(D)
    except TypeError:
      n = 1
    Dbar = np.mean(D)
    invV_0 = 1./self.V_0
    V_n = 1./(invV_0 + n)
    m_n = V_n*(invV_0*self.m_0 + n*Dbar)
    a_n = self.a_0 + n/2.0
    b_n = self.b_0 + 0.5*(np.sum((D-Dbar)**2)+n / (1.0+n*self.V_0)*(self.m_0-Dbar)**2)
    return m_n, V_n, a_n, b_n

  def pred(self, x):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """
    return t_density(2.0*self.a_0, self.m_0, self.b_0*(1.0+self.V_0)/self.a_0, x)

  def evidence(self, D):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """
    m_n, V_n, a_n, b_n = self._post_params(D)
    try:
      n = len(D)
    except:
      n = 1
    return (np.sqrt(np.abs(V_n/self.V_0)) * (self.b_0**self.a_0)/(b_n**a_n) *
            gamma(a_n)/gamma(self.a_0) / (np.pi**(n/2.0)*2.0**(n/2.0)))

  def marginal_var(self, var):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """
    nu_0 = 2*self.a_0
    sigsqr_0 = 2*self.b_0/nu_0
    return scaled_IX_density(nu_0, sigsqr_0, var)

  def marginal_mu(self, mu):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """
    mu_0 = self.m_0
    kappa_0 = 1./self.V_0
    nu_0 = 2*self.a_0
    sigsqr_0 = 2*self.b_0/nu_0
    return t_density(nu_0, mu_0, sigsqr_0/kappa_0, mu)



class NormInvWish(Prior):
  """This class represents TODO.

  *Keyword arguments:*

    - argument1 -- Short description. This argument represents a long description. It can be:
      - Possibility 1: A possibility 1.
      - Possibility 2: A possibility 2.

    - argument2 -- Short description. This argument represents a long description. It can be:
      - Possibility 1: A possibility 1.
      - Possibility 2: A possibility 2.
  """

  def __init__(self, mu_0, kappa_0, Lam_0, nu_0):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """
    self.mu_0 = np.array(mu_0, dtype=float)
    self.kappa_0 = float(kappa_0)
    self.Lam_0 = np.array(Lam_0, dtype=float)
    self.nu_0 = int(nu_0)
    self.d = len(mu_0)
    self.model_dtype = np.dtype([('mu', float, self.d), ('Sig', float, (self.d, self.d))])
    super(NormInvWish, self).__init__()

  def _S(self, D):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """
    Dbar = np.mean(D, axis=0)
    return np.dot((D-Dbar).T, (D-Dbar))

  def sample(self, size=None):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """
    Sig = random_invwish(dof=self.nu_0, invS=self.Lam_0, size=size)
    if size is None:
      ret = np.zeros(1, dtype=self.model_dtype)
      ret['Sig'] = Sig
      ret['mu'] = np.random.multivariate_normal(self.mu_0, Sig/self.kappa_0)
      return ret[0]
    else:
      ret = np.zeros(size, dtype=self.model_dtype)
      ret['Sig'] = Sig
      for r in ret.ravel():
        r['mu'] = np.random.multivariate_normal(self.mu_0, r['Sig']/self.kappa_0)
      return ret

  def like1(self, *args):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """
    if len(args) == 2:
      x, theta = args
      mu = theta['mu']
      Sig = theta['Sig']
    elif len(args) == 3:
      x, mu, Sig = args
    assert x.shape[-1] == self.d
    assert mu.shape[-1] == self.d
    assert Sig.shape[-1] == Sig.shape[-2] == self.d
    norm = np.sqrt((2*np.pi)**self.d * np.linalg.det(Sig))
    einsum = np.einsum("...i,...ij,...j", x-mu, np.linalg.inv(Sig), x-mu)
    return np.exp(-0.5*einsum)/norm

  def __call__(self, *args):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """
    if len(args) == 1:
      mu = args[0]['mu']
      Sig = args[0]['Sig']
    elif len(args) == 2:
      mu, Sig = args
    nu_0, d = self.nu_0, self.d
    Z = (2.0**(nu_0*d/2.0) * gammad(d, nu_0/2.0) *
        (2.0*np.pi/self.kappa_0)**(d/2.0) / np.linalg.det(self.Lam_0)**(nu_0/2.0))
    detSig = np.linalg.det(Sig)
    invSig = np.linalg.inv(Sig)
    einsum = np.einsum("...i,...ij,...j", mu-self.mu_0, invSig, mu-self.mu_0)
    return 1./Z * detSig**(-((nu_0+d)/2.0+1.0)) *
           np.exp(-0.5*np.trace(np.einsum("...ij,...jk->...ik", self.Lam_0, invSig), axis1=-2, axis2=-1) - self.kappa_0/2.0*einsum)

  def _post_params(self, D):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """
    shape = D.shape
    if len(shape) == 2:
      n = shape[0]
      Dbar = np.mean(D, axis=0)
    elif len(shape) == 1:
      n = 1
      Dbar = np.mean(D)
    kappa_n = self.kappa_0 + n
    nu_n = self.nu_0 + n
    mu_n = (self.kappa_0 * self.mu_0 + n * Dbar) / kappa_n
    x = (Dbar-self.mu_0)[:, np.newaxis]
    Lam_n = (self.Lam_0 + self._S(D) + self.kappa_0*n/kappa_n*np.dot(x, x.T))
    return mu_n, kappa_n, Lam_n, nu_n

  def pred(self, x):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """
    return multivariate_t_density(self.nu_0-self.d+1, self.mu_0, self.Lam_0*(self.kappa_0+1)/(self.kappa_0 - self.d + 1), x)

  def evidence(self, D):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """
    shape = D.shape
    if len(shape) == 2:
      n, d = shape
    elif len(shape) == 1:
      n, d = 1, shape[0]
    assert d == self.d
    mu_n, kappa_n, Lam_n, nu_n = self._post_params(D)
    detLam0 = np.linalg.det(self.Lam_0)
    detLamn = np.linalg.det(Lam_n)
    num = gammad(d, nu_n/2.0) * detLam0**(self.nu_0/2.0)
    den = np.pi**(n*d/2.0) * gammad(d, self.nu_0/2.0) * detLamn**(nu_n/2.0)
    return num/den * (self.kappa_0/kappa_n)**(d/2.0)




















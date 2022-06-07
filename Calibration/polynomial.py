from __future__ import division, absolute_import, print_function

import warnings
import numpy as np
import numpy.linalg as la

from numpy.polynomial.polynomial import polyvander, polyvander2d, polyvander3d
from numpy.polynomial import polyutils as pu

# inspired from numpy.polynomial.polynomial.polyfit DV 7/2014
def polyfit2d(x, y, val, deg, rcond=None, full=False, w=None):
    """
    Least-squares fit of a 2d polynomial to 3d data.

    Return the coefficients of a polynomial of degree `deg` that is the
    least squares fit to the data values `val` given at points `x,y`. 
    If `val` is 1-D the returned coefficients will also be 1-D. If `val` 
    is 2-D multiple fits are done, one for each column of `val`, and the 
    resulting coefficients are stored in the corresponding columns of a 2-D 
    return.
    The fitted polynomial(s) are in the form

    .. math::  p(x,y) = \\sum_{i,j} c_{i,j} * x^i * y^j
    
    where `0 <= i <= l`, `0 <= j <= m`, and `l, m = deg` 
    are the given degrees in `x, y`. 

    Since numpy version 1.7.0, polyfit also supports NA. If any of the
    elements of `x`, `y`, 'val' or `w` are NA, then the corresponding 
    rows of the linear least squares problem (see Notes) are set to 0. 
    If `val` is 2-D, then an NA in any row of `val` invalidates that whole row.

    Parameters
    ----------
    x, y : array_like, shape (`M`,)
        x,y-coordinates of the `M` sample (data) points ``(x[i], y[i], val[i])``.
    val : array_like, shape (`M`,) or (`M`, `K`)
        values of to be fitted at the sample points.  Several sets of sample points
        sharing the same x,y-coordinates can be (independently) fit with one
        call to `polyfit2d` by passing in for `val` a 2-D array that contains
        one data set per column.
    deg : list of ints
        List of maximum degrees of the form [x_deg, y_deg]
    rcond : float, optional
        Relative condition number of the fit.  Singular values smaller
        than `rcond`, relative to the largest singular value, will be
        ignored.  The default value is ``len(x)*eps``, where `eps` is the
        relative precision of the platform's float type, about 2e-16 in
        most cases.
    full : bool, optional
        Switch determining the nature of the return value.  When ``False``
        (the default) just the coefficients are returned; when ``True``,
        diagnostic information from the singular value decomposition (used
        to solve the fit's matrix equation) is also returned.
    w : array_like, shape (`M`,), optional
        Weights. If not None, the contribution of each point
        ``(x[i],y[i],val[i])`` to the fit is weighted by `w[i]`. Ideally 
        the weights are chosen so that the errors of the products ``w[i]*val[i]``
        all have the same variance.  The default value is None.


    Returns
    -------
    coef : ndarray, (`order`,) or (`order`, `K`), where
        :math:`order = (deg[0]+1)*(deg([1]+1)`.
        Polynomial coefficients ordered from low to high.  If `val` was 2-D,
        the coefficients in column `k` of `coef` represent the polynomial
        fit to the data in `val`'s `k`-th column.

    [residuals, rank, singular_values, rcond] : present when `full` == True
        Sum of the squared residuals (SSR) of the least-squares fit; the
        effective rank of the scaled Vandermonde matrix; its singular
        values; and the specified value of `rcond`.  For more information,
        see `linalg.lstsq`.

    Raises
    ------
    RankWarning
        Raised if the matrix in the least-squares fit is rank deficient.
        The warning is only raised if `full` == False.  The warnings can
        be turned off by:

        >>> import warnings
        >>> warnings.simplefilter('ignore', RankWarning)

    See Also
    --------
    chebfit, legfit, lagfit, hermfit, hermefit
    polyfit
    polyval, polyval2d : Evaluates a polynomial.
    polyvander, polyvander2d : Vandermonde matrix for powers.
    linalg.lstsq : Computes a least-squares fit from the matrix.
    scipy.interpolate.UnivariateSpline : Computes spline fits.

    Notes
    -----
    The solution is the coefficients of the polynomial `p` that minimizes
    the sum of the weighted squared errors

    .. math :: E = \\sum_j w_j^2 * |val_j - p(x_j,y_j)|^2,

    where the :math:`w_j` are the weights. This problem is solved by
    setting up the (typically) over-determined matrix equation:

    .. math :: V(x,y) * c = w * val,

    where `V` is the weighted pseudo  Vandermonde matrix of `x,y`, `c` are the
    coefficients to be solved for, `w` are the weights, and `val` are the
    observed values.  This equation is then solved using the singular value
    decomposition of `V`.

    If some of the singular values of `V` are so small that they are
    neglected (and `full` == ``False``), a `RankWarning` will be raised.
    This means that the coefficient values may be poorly determined.
    Fitting to a lower order polynomial will usually get rid of the warning
    (but may not be what you want, of course; if you have independent
    reason(s) for choosing the degree which isn't working, you may have to:
    a) reconsider those reasons, and/or b) reconsider the quality of your
    data).  The `rcond` parameter can also be set to a value smaller than
    its default, but the resulting fit may be spurious and have large
    contributions from roundoff error.

    Polynomial fits using double precision tend to "fail" at about
    (polynomial) degree 20. Fits using Chebyshev or Legendre series are
    generally better conditioned, but much can still depend on the
    distribution of the sample points and the smoothness of the data.  If
    the quality of the fit is inadequate, splines may be a good
    alternative.

    Examples
    --------
    >>> from numpy import polynomial as P
    >>> x = np.linspace(-1,1,51) # x "data": [-1, -0.96, ..., 0.96, 1]
    >>> y = x.copy()
    >>> x,y = np.meshgrid(x,y,indexing='ij')
    >>> x = x.ravel()
    >>> y = y.ravel()
    >>> val = x**3*y**3 - x + np.random.randn(len(x)) 
    >>> c, stats = polyfit2d(x,y,val,[3,3],full=True)
    >>> c.shape=(4,4)
    >>> c[3,3],c[1,0] # should be approx 1.,-1.
     (0.5958946702829564, -0.92927870426527537)
    
    Same thing without the added noise

    >>> val = x**3*y**3 - x 
    >>> c, stats = polyfit2d(x,y,val,[3,3],full=True)
    >>> c.shape=(4,4)
    >>> c[3,3],c[1,0]
    (0.99999999999999567, -1.0000000000000002)

    """
    ideg = [int(d) for d in deg]
    order = np.prod(np.array(ideg) + 1)
    x, y = np.array((x, y), copy=0, dtype=float)
    val = np.asarray(val, dtype=float)
    
    # check arguments.
    is_valid = [id == d and id >= 0 for id, d in zip(ideg, deg)]
    if is_valid != [1, 1]:
        raise ValueError("degrees must be non-negative integers")    
    if x.ndim != 1 or y.ndim != 1:
        raise TypeError("expected 1D vector for x and y")
    if x.size == 0 or y.size == 0:
        raise TypeError("expected non-empty vector for x and y")
    if val.ndim < 1 or val.ndim > 2 :
        raise TypeError("expected 1D or 2D array for val")
    if not np.all(len(val) == np.array([len(x), len(y)])):
        raise TypeError("expected x,y and val to have same length")

    # set up the least squares matrices in transposed form
    lhs = polyvander2d(x, y, deg).T
    rhs = val.T
    if w is not None:
        w = np.asarray(w) + 0.0
        if w.ndim != 1:
            raise TypeError("expected 1D vector for w")
        if len(x) != len(w):
            raise TypeError("expected x and w to have same length")
        # apply weights. Don't use inplace operations as they
        # can cause problems with NA.
        lhs = lhs * w
        rhs = rhs * w

    # set rcond
    if rcond is None :
        rcond = len(x)*np.finfo(x.dtype).eps

    # Determine the norms of the design matrix columns.
    if issubclass(lhs.dtype.type, np.complexfloating):
        scl = np.sqrt((np.square(lhs.real) + np.square(lhs.imag)).sum(1))
    else:
        scl = np.sqrt(np.square(lhs).sum(1))
    scl[scl == 0] = 1

    # Solve the least squares problem.
    c, resids, rank, s = la.lstsq(lhs.T/scl, rhs.T, rcond)
    c = (c.T/scl).T
    if val.ndim == 1:   
        c.shape = tuple(np.array(deg) + 1)
    else:
        c.shape = tuple(np.array(deg) + 1) + (val.shape[1],)

    # warn on rank reduction
    if rank != order and not full:
        msg = "The fit may be poorly conditioned"
        warnings.warn(msg, pu.RankWarning)

    if full :
        return c, [resids, rank, s, rcond]
    else :
        return c


# inspired from numpy.polynomial.polynomial.polyfit DV 7/2014
def polyfit3d(x, y, z, val, deg, rcond=None, full=False, w=None):
    """
    Least-squares fit of a 3d polynomial to 3d data.

    Return the coefficients of a polynomial of degree `deg` that is the
    least squares fit to the data values `val` given at points `x,y,z`. 
    If `val` is 1-D the returned coefficients will also be 1-D. If `val` 
    is 2-D multiple fits are done, one for each column of `val`, and the 
    resulting coefficients are stored in the corresponding columns of a 2-D 
    return.
    The fitted polynomial(s) are in the form

    .. math::  p(x,y,z) = \\sum_{i,j,k} c_{i,j,k} * x^i * y^j * z^k
    
    where `0 <= i <= l`, `0 <= j <= m`, and `0 <= k <= n`, and `l, m, n = deg` 
    are the given degrees in `x, y, z`. 

    Since numpy version 1.7.0, polyfit also supports NA. If any of the
    elements of `x`, `y`, 'z', 'val' or `w` are NA, then the corresponding 
    rows of the linear least squares problem (see Notes) are set to 0. 
    If `val` is 2-D, then an NA in any row of `val` invalidates that whole row.

    Parameters
    ----------
    x, y, z : array_like, shape (`M`,)
        x,y,z-coordinates of the `M` sample (data) points ``(x[i], y[i], z[i], val[i])``.
    val : array_like, shape (`M`,) or (`M`, `K`)
        values of to be fitted at the sample points.  Several sets of sample points
        sharing the same x,y,z-coordinates can be (independently) fit with one
        call to `polyfit3d` by passing in for `val` a 2-D array that contains
        one data set per column.
    deg : list of ints
        List of maximum degrees of the form [x_deg, y_deg, z_deg]
    rcond : float, optional
        Relative condition number of the fit.  Singular values smaller
        than `rcond`, relative to the largest singular value, will be
        ignored.  The default value is ``len(x)*eps``, where `eps` is the
        relative precision of the platform's float type, about 2e-16 in
        most cases.
    full : bool, optional
        Switch determining the nature of the return value.  When ``False``
        (the default) just the coefficients are returned; when ``True``,
        diagnostic information from the singular value decomposition (used
        to solve the fit's matrix equation) is also returned.
    w : array_like, shape (`M`,), optional
        Weights. If not None, the contribution of each point
        ``(x[i],y[i],z[i],val[i])`` to the fit is weighted by `w[i]`. Ideally 
        the weights are chosen so that the errors of the products ``w[i]*val[i]``
        all have the same variance.  The default value is None.


    Returns
    -------
    coef : ndarray, (`order`,) or (`order`, `K`), where
        :math:`order = (deg[0]+1)*(deg([1]+1)*(deg[2]+1)`.
        Polynomial coefficients ordered from low to high.  If `val` was 2-D,
        the coefficients in column `k` of `coef` represent the polynomial
        fit to the data in `val`'s `k`-th column.

    [residuals, rank, singular_values, rcond] : present when `full` == True
        Sum of the squared residuals (SSR) of the least-squares fit; the
        effective rank of the scaled Vandermonde matrix; its singular
        values; and the specified value of `rcond`.  For more information,
        see `linalg.lstsq`.

    Raises
    ------
    RankWarning
        Raised if the matrix in the least-squares fit is rank deficient.
        The warning is only raised if `full` == False.  The warnings can
        be turned off by:

        >>> import warnings
        >>> warnings.simplefilter('ignore', RankWarning)

    See Also
    --------
    chebfit, legfit, lagfit, hermfit, hermefit
    polyfit
    polyval, polyval3d : Evaluates a polynomial.
    polyvander, polyvander3d : Vandermonde matrix for powers.
    linalg.lstsq : Computes a least-squares fit from the matrix.
    scipy.interpolate.UnivariateSpline : Computes spline fits.

    Notes
    -----
    The solution is the coefficients of the polynomial `p` that minimizes
    the sum of the weighted squared errors

    .. math :: E = \\sum_j w_j^2 * |val_j - p(x_j,y_j,z_j)|^2,

    where the :math:`w_j` are the weights. This problem is solved by
    setting up the (typically) over-determined matrix equation:

    .. math :: V(x,y,z) * c = w * val,

    where `V` is the weighted pseudo  Vandermonde matrix of `x,y,z`, `c` are the
    coefficients to be solved for, `w` are the weights, and `val` are the
    observed values.  This equation is then solved using the singular value
    decomposition of `V`.

    If some of the singular values of `V` are so small that they are
    neglected (and `full` == ``False``), a `RankWarning` will be raised.
    This means that the coefficient values may be poorly determined.
    Fitting to a lower order polynomial will usually get rid of the warning
    (but may not be what you want, of course; if you have independent
    reason(s) for choosing the degree which isn't working, you may have to:
    a) reconsider those reasons, and/or b) reconsider the quality of your
    data).  The `rcond` parameter can also be set to a value smaller than
    its default, but the resulting fit may be spurious and have large
    contributions from roundoff error.

    Polynomial fits using double precision tend to "fail" at about
    (polynomial) degree 20. Fits using Chebyshev or Legendre series are
    generally better conditioned, but much can still depend on the
    distribution of the sample points and the smoothness of the data.  If
    the quality of the fit is inadequate, splines may be a good
    alternative.

    Examples
    --------
    >>> from numpy import polynomial as P
    >>> x = np.linspace(-1,1,51) # x "data": [-1, -0.96, ..., 0.96, 1]
    >>> y = x.copy()
    >>> z = x.copy()
    >>> x,y,z = np.meshgrid(x,y,z,indexing='ij')
    >>> x = x.ravel()
    >>> y = y.ravel()
    >>> z = z.ravel()
    >>> val = x**3*y**3 - x + z + np.random.randn(len(x)) 
    >>> c, stats = polyfit3d(x,y,z,val,[3,3,3],full=True)
    >>> c.shape=(4,4,4)
    >>> c[3,3,0],c[1,0,0],c[0,0,1] # should be approx 1.,-1.,1.
    (1.0273207730911684, -1.0277688088223322, 1.0128996797552634)
    
    Same thing without the added noise

    >>> y = x**3*y**3 - x + z
    >>> c, stats = polyfit3d(x,y,z,val,[3,3,3],full=True)
    >>> c.shape=(4,4,4)
    >>> c[3,3,0],c[1,0,0],c[0,0,1]
    (1.0000000000000102, -0.99999999999998335, 1.00000000000004)

    """
    ideg = [int(d) for d in deg]
    order = np.prod(np.array(ideg) + 1)
    x, y, z = np.array((x, y, z), copy=0, dtype=float)
    val = np.asarray(val, dtype=float)
    
    # check arguments.
    is_valid = [id == d and id >= 0 for id, d in zip(ideg, deg)]
    if is_valid != [1, 1, 1]:
        raise ValueError("degrees must be non-negative integers")    
    if x.ndim != 1 or y.ndim != 1 or z.ndim != 1:
        raise TypeError("expected 1D vector for x,y and z")
    if x.size == 0 or y.size == 0 or z.size == 0:
        raise TypeError("expected non-empty vector for x,y and z")
    if val.ndim < 1 or val.ndim > 2 :
        raise TypeError("expected 1D or 2D array for val")
    if not np.all(len(val) == np.array([len(x), len(y), len(z)])):
        raise TypeError("expected x,y,z and val to have same length")

    # set up the least squares matrices in transposed form
    lhs = polyvander3d(x, y, z, deg).T
    rhs = val.T
    if w is not None:
        w = np.asarray(w) + 0.0
        if w.ndim != 1:
            raise TypeError("expected 1D vector for w")
        if len(x) != len(w):
            raise TypeError("expected x and w to have same length")
        # apply weights. Don't use inplace operations as they
        # can cause problems with NA.
        lhs = lhs * w
        rhs = rhs * w

    # set rcond
    if rcond is None :
        rcond = len(x)*np.finfo(x.dtype).eps

    # Determine the norms of the design matrix columns.
    if issubclass(lhs.dtype.type, np.complexfloating):
        scl = np.sqrt((np.square(lhs.real) + np.square(lhs.imag)).sum(1))
    else:
        scl = np.sqrt(np.square(lhs).sum(1))
    scl[scl == 0] = 1

    # Solve the least squares problem.
    c, resids, rank, s = la.lstsq(lhs.T/scl, rhs.T, rcond)
    c = (c.T/scl).T
    if val.ndim == 1:   
        c.shape = tuple(np.array(deg) + 1)
    else:
        c.shape = tuple(np.array(deg) + 1) + (val.shape[1],)

    # warn on rank reduction
    if rank != order and not full:
        msg = "The fit may be poorly conditioned"
        warnings.warn(msg, pu.RankWarning)

    if full :
        return c, [resids, rank, s, rcond]
    else :
        return c




# inspired from numpy.polynomial.polynomial.polyfit DV 11/2014
def polyfit2d_wocty(x, y, val, deg, rcond=None, full=False, w=None):
    """
    Least-squares fit of a 2d polynomial (with no constant term) to 3d data.

    Return the coefficients of a polynomial of degree `deg` that is the
    least squares fit to the data values `val` given at points `x,y`. 
    If `val` is 1-D the returned coefficients will also be 1-D. If `val` 
    is 2-D multiple fits are done, one for each column of `val`, and the 
    resulting coefficients are stored in the corresponding columns of a 2-D 
    return.
    The fitted polynomial(s) are in the form

    .. math::  p(x,y) = \\sum_{i,j} c_{i,j} * x^i * y^j 
    
    where `0 <= i <= l`, `0 <= j <= m`, and `l, m = deg` 
    are the given degrees in `x, y`. 

    Since numpy version 1.7.0, polyfit also supports NA. If any of the
    elements of `x`, `y`, 'val' or `w` are NA, then the corresponding 
    rows of the linear least squares problem (see Notes) are set to 0. 
    If `val` is 2-D, then an NA in any row of `val` invalidates that whole row.

    Parameters
    ----------
    x, y : array_like, shape (`M`,)
        x,y-coordinates of the `M` sample (data) points ``(x[i], y[i], val[i])``.
    val : array_like, shape (`M`,) or (`M`, `K`)
        values of to be fitted at the sample points.  Several sets of sample points
        sharing the same x,y-coordinates can be (independently) fit with one
        call to `polyfit2d` by passing in for `val` a 2-D array that contains
        one data set per column.
    deg : list of ints
        List of maximum degrees of the form [x_deg, y_deg]
    rcond : float, optional
        Relative condition number of the fit.  Singular values smaller
        than `rcond`, relative to the largest singular value, will be
        ignored.  The default value is ``len(x)*eps``, where `eps` is the
        relative precision of the platform's float type, about 2e-16 in
        most cases.
    full : bool, optional
        Switch determining the nature of the return value.  When ``False``
        (the default) just the coefficients are returned; when ``True``,
        diagnostic information from the singular value decomposition (used
        to solve the fit's matrix equation) is also returned.
    w : array_like, shape (`M`,), optional
        Weights. If not None, the contribution of each point
        ``(x[i],y[i],val[i])`` to the fit is weighted by `w[i]`. Ideally 
        the weights are chosen so that the errors of the products ``w[i]*val[i]``
        all have the same variance.  The default value is None.


    Returns
    -------
    coef : ndarray, (`order`,) or (`order`, `K`), where
        :math:`order = (deg[0]+1)*(deg([1])`.
        Polynomial coefficients ordered from low to high.  If `val` was 2-D,
        the coefficients in column `k` of `coef` represent the polynomial
        fit to the data in `val`'s `k`-th column.

    [residuals, rank, singular_values, rcond] : present when `full` == True
        Sum of the squared residuals (SSR) of the least-squares fit; the
        effective rank of the scaled Vandermonde matrix; its singular
        values; and the specified value of `rcond`.  For more information,
        see `linalg.lstsq`.

    Raises
    ------
    RankWarning
        Raised if the matrix in the least-squares fit is rank deficient.
        The warning is only raised if `full` == False.  The warnings can
        be turned off by:

        >>> import warnings
        >>> warnings.simplefilter('ignore', RankWarning)

    See Also
    --------
    chebfit, legfit, lagfit, hermfit, hermefit
    polyfit
    polyval, polyval3d : Evaluates a polynomial.
    polyvander, polyvander3d : Vandermonde matrix for powers.
    linalg.lstsq : Computes a least-squares fit from the matrix.
    scipy.interpolate.UnivariateSpline : Computes spline fits.

    Notes
    -----
    The solution is the coefficients of the polynomial `p` that minimizes
    the sum of the weighted squared errors

    .. math :: E = \\sum_j w_j^2 * |val_j - p(x_j,y_j)|^2,

    where the :math:`w_j` are the weights. This problem is solved by
    setting up the (typically) over-determined matrix equation:

    .. math :: V(x,y) * c = w * val,

    where `V` is the weighted pseudo  Vandermonde matrix of `x,y`, `c` are the
    coefficients to be solved for, `w` are the weights, and `val` are the
    observed values.  This equation is then solved using the singular value
    decomposition of `V`.

    If some of the singular values of `V` are so small that they are
    neglected (and `full` == ``False``), a `RankWarning` will be raised.
    This means that the coefficient values may be poorly determined.
    Fitting to a lower order polynomial will usually get rid of the warning
    (but may not be what you want, of course; if you have independent
    reason(s) for choosing the degree which isn't working, you may have to:
    a) reconsider those reasons, and/or b) reconsider the quality of your
    data).  The `rcond` parameter can also be set to a value smaller than
    its default, but the resulting fit may be spurious and have large
    contributions from roundoff error.

    Polynomial fits using double precision tend to "fail" at about
    (polynomial) degree 20. Fits using Chebyshev or Legendre series are
    generally better conditioned, but much can still depend on the
    distribution of the sample points and the smoothness of the data.  If
    the quality of the fit is inadequate, splines may be a good
    alternative.

    Examples
    --------
    >>> from numpy import polynomial as P
    >>> x = np.linspace(-1,1,51) # x "data": [-1, -0.96, ..., 0.96, 1]
    >>> y = x.copy()
    >>> x,y = np.meshgrid(x,y,indexing='ij')
    >>> x = x.ravel()
    >>> y = y.ravel()
    >>> val = x**3*y**3 - x + np.random.randn(len(x)) 
    >>> c, stats = polyfit2d(x,y,val,[3,3],full=True)
    >>> c.shape=(4,4)
    
    Same thing without the added noise

    >>> y = x**3 - x 
    >>> c, stats = polyfit2d(x,y,val,[3,3],full=True)
    >>> c.shape=(4,4)

    """
    ideg = [int(d) for d in deg]
    order = np.prod(np.array(ideg) + [1, 0]) # remove one because we supress the cst
    x, y = np.array((x, y), dtype=float, copy=0)
    val = np.asarray(val, dtype=float)
    
    # check arguments.
    is_valid = [id == d and id >= 0 for id, d in zip(ideg, deg)]
    if is_valid != [1, 1]:
        raise ValueError("degrees must be non-negative integers")    
    if x.ndim != 1 or y.ndim != 1:
        raise TypeError("expected 1D vector for x and y")
    if x.size == 0 or y.size == 0:
        raise TypeError("expected non-empty vector for x and y")
    if val.ndim < 1 or val.ndim > 2 :
        raise TypeError("expected 1D or 2D array for val")
    if not np.all(len(val) == np.array([len(x), len(y)])):
        raise TypeError("expected x,y and val to have same length")

    # set up the least squares matrices in transposed form
    # removing the 1st column of 1 corresponding to the constant term (x^0*y^0)
    lhs = polyvander2d_wocty(x, y, deg).T
    rhs = val.T
    if w is not None:
        w = np.asarray(w) + 0.0
        if w.ndim != 1:
            raise TypeError("expected 1D vector for w")
        if len(x) != len(w):
            raise TypeError("expected x and w to have same length")
        # apply weights. Don't use inplace operations as they
        # can cause problems with NA.
        lhs = lhs * w
        rhs = rhs * w

    # set rcond
    if rcond is None :
        rcond = len(x)*np.finfo(x.dtype).eps

    # Determine the norms of the design matrix columns.
    if issubclass(lhs.dtype.type, np.complexfloating):
        scl = np.sqrt((np.square(lhs.real) + np.square(lhs.imag)).sum(1))
    else:
        scl = np.sqrt(np.square(lhs).sum(1))
    scl[scl == 0] = 1

    # Solve the least squares problem.
    c, resids, rank, s = la.lstsq(lhs.T/scl, rhs.T, rcond)
    c = (c.T/scl).T
    
    # set the constant coeff to zero
    #c = np.hstack([np.zeros(c[:,[0],...].shape), c])

    # reshape c to degx+1, degy (for use by polyval2d) 
    if val.ndim == 1: 
        c.shape = tuple(np.array(deg) + [1,0]) 
    else:
        c.shape = tuple(np.array(deg) + [1,0]) + (val.shape[1],)
    
    # warn on rank reduction
    if rank != order and not full:
        msg = "The fit may be poorly conditioned"
        warnings.warn(msg, pu.RankWarning)

    if full :
        return c, [resids, rank, s, rcond]
    else :
        return c
    

def polyvander2d_wocty(x, y, deg):
    """Pseudo-Vandermonde matrix of given degrees.

    Returns the pseudo-Vandermonde matrix of degrees `deg` and sample
    points `(x, y)`. The pseudo-Vandermonde matrix is defined by

    .. math:: V[..., deg[1]*i + j] = x^i * y^j,

    where `0 <= i <= deg[0]` and `0 <= j <= deg[1]`. The leading indices of
    `V` index the points `(x, y)` and the last index encodes the powers of
    `x` and `y`.

    If ``V = polyvander2d(x, y, [xdeg, ydeg])``, then the columns of `V`
    correspond to the elements of a 2-D coefficient array `c` of shape
    (xdeg + 1, ydeg + 1) in the order

    .. math:: c_{00}, c_{01}, c_{02} ... , c_{10}, c_{11}, c_{12} ...

    and ``np.dot(V, c.flat)`` and ``polyval2d(x, y, c)`` will be the same
    up to roundoff. This equivalence is useful both for least squares
    fitting and for the evaluation of a large number of 2-D polynomials
    of the same degrees and sample points.

    Parameters
    ----------
    x, y : array_like
        Arrays of point coordinates, all of the same shape. The dtypes
        will be converted to either float64 or complex128 depending on
        whether any of the elements are complex. Scalars are converted to
        1-D arrays.
    deg : list of ints
        List of maximum degrees of the form [x_deg, y_deg].

    Returns
    -------
    vander2d : ndarray
        The shape of the returned matrix is ``x.shape + (order,)``, where
        :math:`order = (deg[0]+1)*(deg([1]+1)`.  The dtype will be the same
        as the converted `x` and `y`.

    See Also
    --------
    polyvander, polyvander3d. polyval2d, polyval3d

    """
    ideg = [int(d) for d in deg]
    is_valid = [id == d and id >= 0 for id, d in zip(ideg, deg)]
    if is_valid != [1, 1]:
        raise ValueError("degrees must be non-negative integers")
    degx, degy = ideg
    x, y = np.array((x, y), copy=0) + 0.0

    vx = polyvander(x, degx)
    vy = polyvander(y, degy)[...,1:]
    v = vx[..., None]*vy[..., None,:]
    # einsum bug
    #v = np.einsum("...i,...j->...ij", vx, vy)
    return v.reshape(v.shape[:-2] + (-1,))

    

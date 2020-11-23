import numpy as np
import theano.tensor as tt
import theano
from pymc3.theanof import floatX
import pymc3
from scipy.linalg import helmert

class sStickBreaking(pymc3.distributions.transforms.transform):
    """
    Transforms K - 1 dimensional simplex space (k values in [0,1] and that sum to 1) to a K - 1 vector of real values.
    """

    name = "stickbreaking"
    
    def __init__(self, p=9, scale=1):
        self.p = p
        self.b = 1/p
        self.scale = scale
        self.unscale = 1/scale

    def forward(self, x_):
        x = x_.T
        n = x.shape[0]
        lx = tt.log(x)
        shift = tt.sum(lx, 0, keepdims=True) / n
        y = lx[:-1] - shift
        y = y * self.scale
        y = tt.switch(tt.le(y, 0), -((-y)**self.p), y**self.p)
        return floatX(y.T)

    def forward_val(self, x_):
        x = x_.T
        n = x.shape[0]
        lx = np.log(x)
        shift = np.sum(lx, 0, keepdims=True) / n
        y = lx[:-1] - shift
        y = y * self.scale
        y = np.where(y<0, -((-y)**self.p), y**self.p)
        return floatX(y.T)

    def backward(self, y_):
        y = y_.T
        y = tt.switch(tt.le(y, 0), -((-y)**self.b), y**self.b)
        y = (self.unscale*y)
        y = tt.concatenate([y, -tt.sum(y, 0, keepdims=True)])
        # "softmax" with vector support and no deprication warning:
        e_y = tt.exp(y - tt.max(y, 0, keepdims=True))
        x = e_y / tt.sum(e_y, 0, keepdims=True)
        return floatX(x.T)

    def backward_val(self, y_):
        y = y_.T
        y = np.where(y<0, -((-y)**self.b), y**self.b)
        y = (self.unscale*y)
        y = np.concatenate([y, -np.sum(y, 0, keepdims=True)])
        x = np.exp(y)/np.sum(np.exp(y), 0, keepdims=True)
        return floatX(x.T)
    
    def jacobian_det(self, x_):
        x = x_.T
        x = (self.unscale*x)
        fa = ((1/self.p)-1)*tt.sum(tt.log(tt.abs_(x)), 0, keepdims=True)
        x = tt.switch(tt.le(x, 0), -((-x)**self.b), x**self.b)
        n = x.shape[0]
        sx = tt.sum(x, 0, keepdims=True)
        r = tt.concatenate([x+sx, tt.zeros(sx.shape)])
        # stable according to: http://deeplearning.net/software/theano_versions/0.9.X/NEWS.html
        sr = tt.log(tt.sum(tt.exp(r), 0, keepdims=True))
        d = tt.log(n) + (n*sx) - (n*sr) + fa + ((n-1)*tt.log(self.b)) + self.unscale
        d = tt.switch(tt.isinf(d), 1e99, d) # is tested on the simplex edge where forward maps to inf
        return tt.sum(d, 0).T


class StickBreaking2(pymc3.distributions.transforms.transform):
    """
    Transforms K - 1 dimensional simplex space (k values in [0,1] and that sum to 1) to a K - 1 vector of real values.
    """

    name = "stickbreaking"

    def forward(self, x_):
        x = x_.T
        n = x.shape[0]
        lx = tt.log(x)
        shift = tt.sum(lx, 0, keepdims=True) / n
        y = lx[:-1] - shift
        return floatX(y.T)

    def forward_val(self, x_):
        x = x_.T
        n = x.shape[0]
        lx = np.log(x)
        shift = np.sum(lx, 0, keepdims=True) / n
        y = lx[:-1] - shift
        return floatX(y.T)

    def backward(self, y_):
        y = y_.T
        y = tt.concatenate([y, -tt.sum(y, 0, keepdims=True)])
        # "softmax" with vector support and no deprication warning:
        e_y = tt.exp(y - tt.max(y, 0, keepdims=True))
        x = e_y / tt.sum(e_y, 0, keepdims=True)
        return floatX(x.T)

    def backward_val(self, y_):
        y = y_.T
        y = np.concatenate([y, -np.sum(y, 0, keepdims=True)])
        x = np.exp(y)/np.sum(np.exp(y), 0, keepdims=True)
        return floatX(x.T)

    def jacobian_det(self, y_):
        y = y_.T
        Km1 = y.shape[0]
        sy = tt.sum(y, 0, keepdims=True)
        r = tt.concatenate([y+sy, tt.zeros(sy.shape)])
        # stable according to: http://deeplearning.net/software/theano_versions/0.9.X/NEWS.html
        sr = tt.log(tt.sum(tt.exp(r), 0, keepdims=True))
        d = tt.log(Km1) + (Km1*sy) - (Km1*sr)
        return tt.sum(d, 0).T



def gram_schmidt(vectors, tol=1e-10):
    basis = []
    for v in vectors:
        w = v - sum(np.dot(v,b)*b for b in basis)
        if (w > tol).any():  
            basis.append(w/np.linalg.norm(w))
    return np.array(basis)


class StickBreaking3(pymc3.distributions.transforms.transform):
    """
    Transforms K - 1 dimensional simplex space (k values in [0,1] and that sum to 1) to a K - 1 vector of real values.
    """

    name = "stickbreaking"

    def __init__(self, n):
        A = np.identity(n)
        A[0, :] = np.ones(n)
        A = gram_schmidt(A)
        self.A = A[1:,:]
        self.tA = theano.shared(A[1:,:])

    def forward(self, x_):
        x = x_.T
        lx = tt.log(x)
        y = tt.dot(self.tA, lx)
        return floatX(y.T)

    def forward_val(self, x_):
        x = x_.T
        lx = np.log(x)
        y = np.dot(self.A, lx)
        return floatX(y.T)

    def backward(self, y_):
        y = tt.dot(y_, self.tA)
        # "softmax" with vector support and no deprication warning:
        e_y = tt.exp(y - tt.max(y, 0, keepdims=True))
        x = e_y / tt.sum(e_y, 0, keepdims=True)
        return floatX(x.T)

    def backward_val(self, y_):
        y = np.dot(y_, self.A)
        x = np.exp(y)/np.sum(np.exp(y), 0, keepdims=True)
        return floatX(x.T)

    def jacobian_det(self, y_):
        y = y_.T
        Km1 = y.shape[0]
        sy = tt.sum(y, 0, keepdims=True)
        r = tt.concatenate([y+sy, tt.zeros(sy.shape)])
        # stable according to: http://deeplearning.net/software/theano_versions/0.9.X/NEWS.html
        sr = tt.log(tt.sum(tt.exp(r), 0, keepdims=True))
        d = tt.log(Km1) + (Km1*sy) - (Km1*sr)
        return tt.sum(d, 0).T




class StickBreaking4(pymc3.distributions.transforms.transform):
    """
    Transforms K - 1 dimensional simplex space (k values in [0,1] and that sum to 1) to a K - 1 vector of real values.
    """

    name = "stickbreaking"

    def __init__(self, n):
        self.A = helmert(n, full=False)
        self.tA = theano.shared(self.A)

    def forward(self, x_):
        x = x_.T
        lx = tt.log(x)
        y = tt.dot(self.tA, lx)
        return floatX(y.T)

    def forward_val(self, x_):
        x = x_.T
        lx = np.log(x)
        y = np.dot(self.A, lx)
        return floatX(y.T)

    def backward(self, y_):
        y = tt.dot(y_, self.tA)
        # "softmax" with vector support and no deprication warning:
        e_y = tt.exp(y - tt.max(y, 0, keepdims=True))
        x = e_y / tt.sum(e_y, 0, keepdims=True)
        return floatX(x.T)

    def backward_val(self, y_):
        y = np.dot(y_, self.A)
        x = np.exp(y)/np.sum(np.exp(y), 0, keepdims=True)
        return floatX(x.T)

    def jacobian_det(self, y_):
        y = y_.T
        Km1 = y.shape[0]
        sy = tt.sum(y, 0, keepdims=True)
        r = tt.concatenate([y+sy, tt.zeros(sy.shape)])
        # stable according to: http://deeplearning.net/software/theano_versions/0.9.X/NEWS.html
        sr = tt.log(tt.sum(tt.exp(r), 0, keepdims=True))
        d = tt.log(Km1) + (Km1*sy) - (Km1*sr)
        return tt.sum(d, 0).T


class StickBreaking5(pymc3.distributions.transforms.transform):
    """
    Transforms K - 1 dimensional simplex space (k values in [0,1] and that sum to 1) to a K - 1 vector of real values.
    """

    name = "stickbreaking"

    def __init__(self, n):
        self.A = helmert(n, full=False)
        self.tA = theano.shared(self.A)

    def forward(self, x_):
        x = x_.T
        lx = tt.log(x)
        y = tt.dot(self.tA, lx)
        return floatX(y.T)

    def forward_val(self, x_):
        x = x_.T
        lx = np.log(x)
        y = np.dot(self.A, lx)
        return floatX(y.T)

    def backward(self, y_):
        y = tt.dot(y_, self.tA)
        # "softmax" with vector support and no deprication warning:
        e_y = tt.exp(y - tt.max(y, 0, keepdims=True))
        x = e_y / tt.sum(e_y, 0, keepdims=True)
        return floatX(x.T)

    def backward_val(self, y_):
        y = np.dot(y_, self.A)
        x = np.exp(y)/np.sum(np.exp(y), 0, keepdims=True)
        return floatX(x.T)

    def jacobian_det(self, y_):
        y = tt.dot(y_, self.tA)
        Km1 = y.shape[0]
        sy = tt.sum(y, 0, keepdims=True)
        r = tt.concatenate([y+sy, tt.zeros(sy.shape)])
        # stable according to: http://deeplearning.net/software/theano_versions/0.9.X/NEWS.html
        sr = tt.log(tt.sum(tt.exp(r), 0, keepdims=True))
        d = tt.log(Km1) + (Km1*sy) - (Km1*sr)
        return tt.sum(d, 0).T


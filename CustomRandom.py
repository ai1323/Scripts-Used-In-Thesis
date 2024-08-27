import numpy as np

class CustomRandom:
    def __init__(self, seed):
        self.state = seed

    def next(self):
        self.state = (self.state * 6364136223846793005 + 1) % 9223372036854775808
        return (self.state / 9223372036854775808) * 0.5773502691896258 - 0.2886751345948129

    def nextVector(self, Hdim):
        vec = np.empty(Hdim, dtype=np.complex128)
        for i in range(Hdim):
            realPart = self.next()
            vec[i] = complex(realPart, 0.0)
        return vec

    def customRandomSeed(self, N, n_skip=1000):
        vec = np.empty(N, dtype=np.int64)
        vec[0] = self.state
        for _ in range(n_skip):
            self.next()
        for i in range(1, N):
            vec[i] = self.state
            self.next()
        return vec.tolist()
    

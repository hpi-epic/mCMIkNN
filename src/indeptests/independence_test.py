from abc import ABC, abstractmethod

class IndependenceTest(ABC):

    @abstractmethod
    def test_params(self):
        pass

    @abstractmethod
    def compute_pval(self, x, y, z):
        pass

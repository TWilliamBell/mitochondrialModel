class namespace: ## For a library solution, it would be possible to use argparse
    ## but my needs are so minimal that I thought it would be simpler to not add
    ## a dependency.
    def __init__(self, **kwargs):
        self.__dict__.update(kwargs)
        self.n_len = len(kwargs.keys())
        self.keys = list(kwargs.keys())
    def len(self):
        return self.n_len

class negativeConcentration(Exception):
    def __init__(self):
        self.message = "Negative concentration appeared."

class outOfOrder(Exception):
    def __init__(self):
        self.message = "Array must be monotonically decreasing."

class timeError(Exception):
    def __init__(self):
        self.message = "DE solver timed out."

class didNotConverge(Exception):
    def __init__(self):
        self.message = "Root-finding method did not succeed."
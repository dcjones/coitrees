class Interval:
    """Interval class for Coitree."""

    def __init__(self, low, high):
        """Initialize interval."""
        self.low = low
        self.high = high


class CoitreeNode:
    """Coitree node class."""

    def __init__(self, interval, data):
        """Initialize node."""
        self.key = interval
        self.data = data


class Coitree:
    """Coitree class."""

    def __init__(self):
        """Initialize Coitree."""
        raise NotImplementedError

    def insert(self, node):
        raise NotImplementedError

    def delete(self, key):
        raise NotImplementedError

    def search(self, key):
        raise NotImplementedError

    def build(self, node_list):
        raise NotImplementedError

    def find_overlap(self, key):
        raise NotImplementedError

    def __repr__(self) -> str:
        raise NotImplementedError

    __str__ = __repr__

# identity.py
# how to identify a function

from .. import IrisCommand
from .. import state_types as t
from .. import iris_objects

class Identity(IrisCommand):
    title = "identity function {x}"
    examples = []
    argument_types = {"x":t.Int("What value?")}
    def command(self, x):
        return x
    def explanation(self, result):
        return result
_Identity = Identity()


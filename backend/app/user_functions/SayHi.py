from iris import state_types as t
from iris import IrisCommand

##### EXAMPLE FUNCTION #######
class SayHi(IrisCommand):
	# what this function will be called/ display name
    title = "Hi"

    # templates of what text can be inputted to call this function
    examples = ["Hello ", "Bonjour ", "Hola", "Sup"]

    # types of inputs want to implement
    argument_types = {}

    # core logic of the command
    def command(self):
        import random
        n0 = random.randint(0,100)
        return n0

    # wrap the output of a command to display to user
    # by default this will be an identity function
    # each element of the list defines a separate chat bubble
    def explanation(self, result):
        return "Your lucky number is "+str(result)

# This line loads the function into the app on startup
_SayHi = SayHi()


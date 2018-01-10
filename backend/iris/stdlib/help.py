# Help.py
# Displays commands of how to use iris and how to enter validation mode


from .. import IrisCommand
from .. import state_types as t
from .. import state_machine as sm
from .. import util as util
from .. import iris_objects



class Help(IrisCommand):
    title = "help"
    examples = ["help", "Help"]
    def command(self):
        return ['Do you have a question you want to ask?', ' Type your question in the input box below', 'If you want to exit a task, type \"stop\"']
    def explanation(self, result):
        return result
_Help = Help()



class Validation(IrisCommand):
	title = "Is this useful?"
	argument_types = {"bool_val":t.YesNo("Is this useful?", yes=True, no=False)}

	def command(self, bool_val):
		#DO SOMESTORAGE
		if bool_val:
			return "Yes"
		else:
			return "No"

	def explanation(self, result):
		return "Exitting validation mode. here's the validation stored"

Validation = Validation()
	# def commands(
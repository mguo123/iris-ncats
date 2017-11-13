# Help.py
# Displays commands of how to use iris


from iris import state_types as t
from iris import IrisCommand

from iris import iris_objects
class Help(IrisCommand):
    title = "help"
    examples = ["help", "Help"]
    def command(self):
        return ['Do you have a question you want to ask?', ' Type your question in the input box below', 'If you want to exit a task, type \"quit\"']
    def explanation(self, result):
        return result
_Help = Help()


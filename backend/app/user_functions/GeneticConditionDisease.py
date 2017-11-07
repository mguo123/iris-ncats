from iris import state_types as t
from iris import IrisCommand

class GeneticConditionDisease(IrisCommand):
    # what iris will call the command + how it will appear in a hint
    title = "how does {condition} protects against {disease}?"
    
    # give an example for iris to recognize the command
    examples = ["find out how {condition} protects against {disease}", "protective mechanism of condition", "protective mechanism of condition against disease", "protective condition", "protective mechanism of {condition} against {disease}", "how does {condition} protect against {disease}"]
    # type annotations for each command argument, to help Iris collect missing values from a user
    argument_types = {"condition":t.String("Just for confirmation: What is the condition you want to analyze?"), "disease":t.String("What is the disease you want to analyze?")}
    
    # core logic of the command
    def command(self, condition, disease):
        import numpy
        return numpy.random.randint(100)
        
    # wrap the output of a command to display to user
    # by default this will be an identity function
    # each element of the list defines a separate chat bubble
    def explanation(self, result):

        return ["Enter the Department of Mysteries to find the answer. Address: Fishbowl, first desk on the right. Here are your magic numbers", result]


_GeneticConditionDisease = GeneticConditionDisease()
from iris import state_types as t
from iris import IrisCommand

class GeneticConditionDisease(IrisCommand):
    # what iris will call the command + how it will appear in a hint
    title = "find out how {genetic_condition} protects against {disease}"
    
    # give an example for iris to recognize the command
    examples = ["find out how {genetic_condition} protects against {disease}", "protective mechanism of genetic_condition", "protective mechanism of genetic_condition against disease", "protective genetic_condition", "protective mechanism of {genetic_condition} against {disease}", "how does {genetic_condition} protect against {disease}"]
    # type annotations for each command argument, to help Iris collect missing values from a user
    argument_types = {"genetic_condition":t.String("What is the genetic_condition you want to analyze?"), "disease":t.String("What is the disease you want to analyze?")}
    
    # core logic of the command
    def command(self, genetic_condition, disease):
        import numpy
        return numpy.random.randint(100)
        
    # wrap the output of a command to display to user
    # by default this will be an identity function
    # each element of the list defines a separate chat bubble
    def explanation(self, result):

        return ["Enter the Department of Mysteries to find the answer. Address: Fishbowl, first desk on the right. Here are your magic numbers", result]


_GeneticConditionDisease = GeneticConditionDisease()
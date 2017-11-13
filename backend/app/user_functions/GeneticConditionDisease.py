from iris import state_types as t
from iris import IrisCommand
from app.user_functions.Q1_main import Q1_query


class GeneticConditionDisease(IrisCommand):
    # what iris will call the command + how it will appear in a hint
    # title = "how does {condition} protects against {condition}?"
    title = "what protects against {condition}"
    # give an example for iris to recognize the command
    examples = ["what protects against {condition}", "what genetic conditions might offer protection against {condition} and why", "protective mechanism of condition against disease", "protective condition", "protection against disease"]
    # type annotations for each command argument, to help Iris collect missing values from a user
    argument_types = {"condition":t.String("Just for confirmation: What is the condition do you want to analyze?")} 
                        # ,"genetic_disease":t.String("What is the genetic disease do you think it might link to? If unknown, type none")}
    
    # core logic of the command
    def command(self, condition):
        # import numpy as np
        return Q1_query(condition)

    # wrap the output of a command to display to user
    # by default this will be an identity function
    # each element of the list defines a separate chat bubble
    def explanation(self, result):

        return ["Enter the Department of Mysteries to find the answer. Address: Fishbowl, first desk on the right. Here are your magic numbers", result]


_GeneticConditionDisease = GeneticConditionDisease()
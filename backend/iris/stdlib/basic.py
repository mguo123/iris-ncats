# basic.py
# Basic numerical, logic, array-based, and arithmetic functions
from .. import IrisCommand
from .. import state_types as t
from .. import state_machine as sm
from .. import util as util
from .. import iris_objects
##############################################################################
############################### Datatype Functions #############################
##############################################################################

class GenerateNumber(IrisCommand):
    title = "generate a random number"
    examples = [ "generate number" ]
    def command(self):
        import random
        return random.randint(0,100)

generateNumber = GenerateNumber()

class GenerateFloat(IrisCommand):
    title = "generate a random float"
    examples = [ "generate float" ]
    def command(self):
        print("I am called by", self.caller.context)
        import random
        return float(random.randint(0,100))

generateFloat = GenerateFloat()

class PrintValue(IrisCommand):
    title = "print {value}"
    examples = [ "display {value}", "{value}", "show me {value}"]
    help_text = [
        "This command will display the underlying data for environment variables."
    ]
    def command(self, value : t.EnvVar()):
        # if isinstance(value, iris_objects.IrisDataframe):
        #     return value.to_matrix()
        return value

printValue = PrintValue()

class ToNumber(IrisCommand):
    title = "{x} to number"
    examples = []
    argument_types = {"x":t.String("What to convert to a number?")}
    def command(self, x):
        try:
            out =  float(x)
        except:
            out = None
        return out
    def explanation(self, result):
        return result
_ToNumber = ToNumber()

##############################################################################
############################### STRING Functions #############################
##############################################################################

class ToString(IrisCommand):
    title = "{x} to string"
    examples = []
    argument_types = {"x":t.Float("What float to convert to string?")}
    def command(self, x):
        return str(x)
    def explanation(self, result):
        return result
_ToString = ToString()

class CountWords(IrisCommand):
    title = "count words in {doc}"
    examples = ["count words {doc}", "word count {doc}"]
    argument_types = {"doc":t.String("What is the string you want to anaylze?")}
    def command(self, doc):
        return len(doc.split())
    def explanation(self, result):
        return "{} words in doc".format(result)
_CountWords = CountWords()

class SplitStringCommas(IrisCommand):
    title = "split string {str} on commas"
    examples = ["splitting string"]
    argument_types = {"str":t.String("What string do you want to split on commas?")}
    def command(self, str):
        return str.split(",")
    def explanation(self, result):
        return [result]
_SplitStringCommas = SplitStringCommas()



##############################################################################
############################### Logic Functions #############################
##############################################################################
class ReturnFalse(IrisCommand):
    title = "return false"
    examples = []
    argument_types = {}
    def command(self, ):
        return False
    def explanation(self, result):
        return result
_ReturnFalse = ReturnFalse()

class ReturnTrue(IrisCommand):
    title = "return true"
    examples = []
    argument_types = {}
    def command(self, ):
        return True
    def explanation(self, result):
        return result
_ReturnTrue = ReturnTrue()

class Equals(IrisCommand):
    title = "{x} equals {y}"
    examples = ["eq {x} {y}", "equal {x} {y}"]
    argument_types = {"x":t.Int("What is the first number?"),"y":t.Int("What is the second number?")}
    def command(self, x, y):
        return x == y
    def explanation(self, result):
        return result
_Equals = Equals()


class LessThan(IrisCommand):
    title = "{x} less than {y}"
    examples = ["{x} < {y}", "{x} less {y}"]
    argument_types = {
        "x": t.Int("Give an integer value for x:"),
        "y": t.Int("Give an integer value for y:")
    }
    def command(self, x, y):
        return x < y

lessThan = LessThan()

class GreaterThan(IrisCommand):
    title = "{x} greater than {y}"
    examples = ["{x} > {y}", "{x} greater {y}"]
    argument_types = {
        "x": t.Int("Give an integer value for x:"),
        "y": t.Int("Give an integer value for y:")
    }
    def command(self, x, y):
        return x > y

greaterThan = GreaterThan()


class CountCharacters(IrisCommand):
    title = "count characters in {string}"
    examples = ["{string} length", "char in {string}"]
    argument_types = {"string":t.String("What string to count characters for?")}
    def command(self, string):
        return len(string)
    def explanation(self, result):
        return result
_CountCharacters = CountCharacters()

##############################################################################
############################### Logic Functions #############################
##############################################################################

class GetArrayLength(IrisCommand):
    title = "get length of array {arr}"
    examples = ["array {arr} length"]
    def command(self, arr: t.Array("What array to get length of?")):
        return arr.shape[0]

getArrayLength = GetArrayLength()

class GenerateArray(IrisCommand):
    title = "generate a random array of {n} numbers"
    examples = ["generate numpy array of size {n}"]
    argument_types = {"n":t.Int("Please enter size of array:")}
    def command(self, n):
        import numpy
        return numpy.random.randint(100, size=n)
    def explanation(self, result):
        return ["Here are the numbers", result]
_GenerateArray = GenerateArray()

class TakeFirst(IrisCommand):
    title = "first {l}"
    examples = ["first element of {list}"]
    argument_types = {"l":t.List("What list?")}
    def command(self, l):
        return l[0]
    def explanation(self, result):
        return result
_TakeFirst = TakeFirst()

class AverageArray(IrisCommand):
    title = "average {array}"
    examples = ["mean {array}"]
    argument_types = {"array":t.Array("What array do you want to average?")}
    def command(self, array):
        import numpy as np
        return np.average(array)
    def explanation(self, result):
        return "The average is {}".format(result)
_AverageArray = AverageArray()



##############################################################################
############################### Arithmetic Functions #############################
##############################################################################
class AddTwo(IrisCommand):
    title = "add two to {x}"
    examples = []
    argument_types = {"x":t.Int("Add two?")}
    def command(self, x):
        return x + 2
    def explanation(self, result):
        return result
_AddTwo = AddTwo()

class AddThree(IrisCommand):
    title = "add three numbers {x} {y} {z}"
    examples = ["add three"]
    argument_types = {"x":t.Int("first int"),"y":t.Int("second int"),"z":t.Int("third int")}
    def command(self, x, y, z):
        return x + y + z
    def explanation(self, result):
        return result
_AddThree = AddThree()

class AddTwoNumbers(IrisCommand):
    title = "add two numbers: {x} and {y}"
    examples = [ "add {x} and {y}",
                 "add {x} {y}",
                 "can you add {x} and {y}" ]
    argument_types = {
        "x": t.Float("Please enter a number for x:"),
        "y": t.Float("Please enter a number for y:")
    }
    help_text = [
        "This command performs addition on two numbers, e.g., 'add 3 and 2' will return 5"
    ]
    def command(self, x, y):
        return x + y

addTwoNumbers = AddTwoNumbers()



class DivideTwoNumbers(IrisCommand):
    title = "divide two numbers: {x} and {y}"
    examples = [ "divide {x} and {y}",
                 "divide {x} {y}" ]
    argument_types = {
        "x": t.Float("Please enter a number for x:"),
        "y": t.Float("Please enter a number for y:")
    }
    help_text = [
        "This command performs addition on two numbers, e.g., 'add 3 and 2' will return 5"
    ]
    def command(self, x, y):
        return x / y
    def explanation(self, val):
        return "Dividing those, I get {}".format(round(val, 6))

divideTwoNumbers = DivideTwoNumbers()


class NegateNumber(IrisCommand):
    title = "negate {x}"
    examples = ["make {x} negative"]
    argument_types = {"x":t.Float("What number to negate?")}
    def command(self, x):
        return x * -1
    def explanation(self, result):
        return result
_NegateNumber = NegateNumber()



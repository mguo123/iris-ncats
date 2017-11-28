'''

Margaret Guo
11/14/17
Create potential modules that NCATS may use: 

Examples include;
How do defects in {xyz pathways} affect {disease phenotype}?
i.e. how do defects in aldehyde processing pathways affect Fanconi Anemia phenotype
What are the risk factors (environmental, genetic, other) that can lead to disease and why?
What are the chances that a {xyz person of given demographics} get {really common disease}
How does {risk factor} affect likelihood to get {disease} [note: based off ehrToPhenotypes]
What are the best techniques I should use to study {biological question}?
Why should people {perform this protective mechanism} (i.e. exercise, drink milk, etc.)?
Basic API calling tasks such as:
Give me the __ ID for this {term}?
What are the pathways associated with these {list of genes}?
Search for databases/best resources
Which database can I use to study {xyz question}?
Which R package can I use to study {xyz question? [biolink]
 
 
'''
from iris import state_types as t
from iris import IrisCommand


from iris import state_machine as sm
from iris import util as util
from iris import iris_objects

class PathwayToPhenotype(IrisCommand):
    title = "How do defects in pathway affect phenotype?"
    examples = ["How do defects in {pathway} affect {phenotype}?", "Why are defs in {pathway} deleterious?", "How is {pathway} involved in the pathogenesis of {phenotype}"]
    argument_types = {"pathway":t.String("What pathway are you interested in?"), "phenotype":t.String("What diseases do you want to analyze?")}
    
    def command(self, pathway, phenotype):
        string_answer = "Currently not implemented, but I also would like to know how defects in %s affect %s" % (pathway, phenotype)
        
        return string_answer

    def explanation(self, result):
        return result

PathwayToPhenotype = PathwayToPhenotype()



class RiskFactorsToDisease(IrisCommand):
    title = "What are the risks (environmental, genetic, etc.) that can lead to disease and why?"
    examples = ["Why is this risk bad for {disease}", "Why are these risks bad for {disease}", "Why should I not do {risks}", "Why should I stopping doing {risks}", "What are the {risks} (environmental, genetic, etc.) that can lead to {disease} and why?"]
    argument_types = {"risks":t.List("What risks do you want to know more about? Enter them in separated by commas"),
                    "disease": t.String("What disease are you analyzing?")}

    def command(self, risks, disease):
        return "Currently not implemented, but that sure is an important question."
        
    def explanation(self, result):
        return result

RiskFactorsToDisease = RiskFactorsToDisease()

class ChanceDisease(IrisCommand):
    title = "What are the chances that a person with conditions gets disease?"

    examples = ["What is the epidemiological risk {conditions} for {disease} and what are their effects on prognosis?", "What is my likelihood of getting {disease}", "How does {conditions} change my chances of getting {disease}?"]
    argument_types = {"conditions":t.List("What conditions (i.e. age, gender, comorbidities, etc.) do you want to factor into this calculation? Please enter them in separated by commas."),
                    "disease": t.String("What disease are you analyzing?")}

    def command(self, conditions, disease):
        return "Currently not implemented:( But we'll get on it as soon as possible."
        
    def explanation(self, result):
        return result    

ChanceDisease = ChanceDisease()


class ResearchQuestion(IrisCommand):
    title = "How can I study {question}?"
    examples = ["What are the best techniques I should use to study {question}?", "Which database can I use to study {question}", "Which R package can I use to study {question}"]
    argument_types = {"question":t.String("What is your research inquiry (clinical, biological, or otherwise)?")}

    def command(self, question):
        return "That's a good question! Unfortunately we aren't magical enough to conjure up the anwer. But we'll definitely get back to you on that."

    def explanation(self, result):
        return result


ResearchQuestion = ResearchQuestion()



'''
Commands to perform basic queries in SmartAPI
'''


from iris import state_types as t
from iris import IrisCommand

from iris import state_machine as sm
from iris import util as util
from iris import iris_objects
from app.user_functions.ncats.tools.SmartAPI import SmartAPI

class ListAllKnowledgeSources(IrisCommand):
    title = "What are all knowledge sources in SmartAPI?"
    examples = ["What are all knowledge sources available?",
                "What available data sources are there?",
                "What resources are available?"]


    def command(self):
        s = SmartAPI.SmartAPI()
        result = s.search_all('*')
        return result

    def explanation(self, result):
        return result


ListAllKnowledgeSources = ListAllKnowledgeSources()

class SearchKnowledgeSourceTitles(IrisCommand):
    title = "What knowledge source titles include {query}?"
    examples = ["What knowledge sources titles contain {query}?"]

    argument_types = {"query": t.String("What is the search term?")}

    def command(self, query):
        s = SmartAPI.SmartAPI()
        result = s.search_titles(query)
        return result # returns list 

    def explanation(self, result):
        if len(result)> 0:
            return result
        else:
            return 'No source titles found'

SearchKnowledgeSourceTitles = SearchKnowledgeSourceTitles()


class SearchKnowledgeSourceFull(IrisCommand):
    title = "What knowledge sources include information about {query}?"
    examples = ["What knowledge sources discuss {query}?",
                "What sources in SmartAPI talk about {query}"]

    argument_types = {"query": t.String("What is the search term?")}

    def command(self, query):
        s = SmartAPI.SmartAPI()
        result = s.search_all(query)
        return result

    def explanation(self, result):
        if len(result)> 0:
            return result
        else:
            return 'No knowledge sources found'

SearchKnowledgeSourceFull = SearchKnowledgeSourceFull()
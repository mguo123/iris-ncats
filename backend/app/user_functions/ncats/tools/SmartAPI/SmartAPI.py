import urllib
import requests
import sys
import pprint


class SmartAPI:
    API_BASE_URL = 'http://smart-api.info/api/query'
    TIMEOUT_SEC = 120

    @staticmethod
    def query(search_term, field=None):
        if field is None:
            url = SmartAPI.API_BASE_URL + '/' + '?' + 'q=' + search_term + '&size=30'
        else:
            search = field + ":" + search_term
            url = SmartAPI.API_BASE_URL + '/' + '?' + 'q=' + search + '&size=30'

        try:
            res = requests.get(url,
                               timeout=SmartAPI.TIMEOUT_SEC)
        except requests.exceptions.Timeout:
            print(url, file=sys.stderr)
            print('Timeout in QueryChEMBL for URL: ' + url, file=sys.stderr)
            return None
        status_code = res.status_code
        if status_code != 200:
            print(url, file=sys.stderr)
            print('Status code ' + str(status_code) + ' for url: ' + url, file=sys.stderr)
            return None
        return res.json()

    @staticmethod
    def search_titles(search_term):
        result = SmartAPI.query(search_term, field='info.title')
        titles = []
        for r in result:
            if r == 'hits':
                #for h in result[r]:
                #     print(h)
                #print(result[r]['info']['title'])
                for h in result[r]:
                    titles.append(h['info']['title'])
        return titles

    @staticmethod
    def search_all(search_term):
        result = SmartAPI.query(search_term)
        titles = []
        for r in result:
            if r == 'hits':
                #for h in result[r]:
                #        print(h)
                    # print(result[r]['info']['title'])
                for h in result[r]:
                    titles.append(h['info']['title'])

        return titles


    #@staticmethod
    #def

if __name__ == '__main__':
    s = SmartAPI
    print(s.search_titles('drug'))
    print(s.search_all('drug'))
    print(s.search_all('*'))









# query parameters paths.parameters:drug  - doesn't work
# query title info.title:gene
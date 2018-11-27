#import urllib2
from six.moves.urllib.request import urlopen,Request
from six.moves.urllib.error import HTTPError
import gzip
import io
import logging
import json
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

class URLResource:

    def __init__(self, url):
        self._url = url
        self._lines = None

    def get_lines(self):
        if self._lines is not None:
            return self._lines
        else:
            try:
                request = Request(self._url)
                request.add_header('User-Agent', 'Mozilla/5.0')
                response = urlopen(request)
            except HTTPError as error:
                logger.error(str(error) + "\n\t " + self._url)
                return []

            raw_content = response.read()
            content = raw_content
            response_info =  response.info()

            if response_info.get('Content-Type') == 'application/json' :
                json_obj = json.loads(content.decode('utf-8'))
                lines = json_obj
            else:

                if response_info.get('Content-Type') == 'application/x-gzip' or \
                    response_info.get('Content-Encoding') == 'gzip' or \
                    self._url.endswith('.gz'):
                    content = gzip.GzipFile(fileobj=io.BytesIO(raw_content)).read()

                content_decoded = content.decode('utf-8')
                #to get lines, split by \n
                lines = content_decoded.split("\n")

            self._lines = lines
        return self._lines

if __name__ == '__main__':
    url_resource = URLResource('http://geneontology.org/ontology/go.obo')
    lines = url_resource.get_lines()
    print("%s\n\t%d lines" % (url_resource._url, len(lines)))
    url_resource = URLResource('ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz')
    lines = url_resource.get_lines()
    print("%s\n\t%d lines" % (url_resource._url, len(lines)))

    url_resource = URLResource('http://www.geneontology.org/gene-associations/goa_human.gaf.json')
    lines = url_resource.get_lines()
    print("%s\n\t%s" % (url_resource._url, lines))

    url_resource = URLResource('http://www.brenda-enzymes.info/ontology/tissue/tree/update/update_files/BrendaTissueOBO')
    lines = url_resource.get_lines()
    print("%s\n\t%s" % (url_resource._url, lines))

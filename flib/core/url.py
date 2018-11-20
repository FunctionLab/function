#import urllib2
from six.moves.urllib.request import urlopen
from six.moves.urllib.error import HTTPError
import gzip
import io
import logging
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
                response = urlopen(self._url)
            except HTTPError as error:
                logger.error(str(error) + "\n\t " + self._url)
                return []

            raw_content = response.read()
            content = raw_content
            response_info =  response.info()
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
    url_resource = URLResource('http://www.geneontology.org/gene-associations/gene_association.human.gz')
    lines = url_resource.get_lines()
    print("%d lines" % len(lines))
    url_resource = URLResource('ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz')
    lines = url_resource.get_lines()
    print("%d lines" % len(lines))




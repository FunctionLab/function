import urllib2
import gzip
import io


class URLResource:

    def __init__(self, url):
        self._url = url
        self._file = None

    def get_file(self):
        if self._file is not None:
            return self._file
        else:
            req = urllib2.Request(self._url, headers={'User-Agent':'Browser'})
            url_file = urllib2.urlopen(req, timeout=10)
            if url_file.info().get('Content-Encoding') == 'gzip' or \
                    url_file.info().get('Content-Type') == 'application/x-gzip':
                url_file = gzip.GzipFile(fileobj=io.BytesIO(url_file.read()))
            self._file = url_file

        return self._file

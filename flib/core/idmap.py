import sys
import logging

logging.basicConfig()
logger = logging.getLogger(__name__)
logger.setLevel(logging.WARNING)


class IDMap:
    """
    Pass the filename of the key_value pair file.
    """

    def __init__(self, filename = None, key_map = None):
        self._key_val = {}
        idfile = list
        if filename is not None:
            idfile = open(filename)

            for line in idfile:
                toks = line.strip().upper().split('\t')
                if len(toks) < 2 or toks[0] == '':
                    continue

                self._key_val[toks[0]] = tuple(toks[1:])
        elif key_map:
            self._key_val = key_map

    def keys(self):
        if self._key_val is None:
            return []
        else:
            return self._key_val.keys()

    """
    Returns None if no file was loaded or if the key does not exist.
    """

    def get(self, id):
        if self._key_val is None:
            logger.info('idmap::get called with no mapping file loaded')
            return None
        else:
            try:
                return self._key_val[id]
            except KeyError:
                logger.warning('No match for %s', id)
                return None

    def __getitem__(self,arg):
        return self.get(arg)

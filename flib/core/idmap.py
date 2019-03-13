import logging

logging.basicConfig()
logger = logging.getLogger(__name__)
logger.setLevel(logging.WARNING)


class IDMap:

    def __init__(self, filename=None, key_map=None):
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
            return list(self._key_val.keys())

    def get(self, id):
        """Returns emtpy list if not loaded or if the key does not exist."""

        if self._key_val is None:
            logger.info('idmap::get called with no mapping file loaded')
            return []
        elif id not in self._key_val:
            logger.warning('No match for %s', id)
            return []
        else:
            return self._key_val[id]

    def __getitem__(self, arg):
        return self.get(arg)

if __name__ == '__main__':
    key_map = {}
    key_map['trnQ'] ={'6775095', '8923208'}
    key_map['MIR10527']={'113218505'}

    idMap = IDMap(key_map=key_map)
    print((list(idMap.keys())))
    print((idMap.get('trnQ')))
    print((idMap.get('xxx')))

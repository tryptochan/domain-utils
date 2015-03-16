class ParseError(Exception) :
    """Super class for exception classes used in Parsers package.
    """

    def __init__(self, file='', msg=''):
        self.file = str(file)
        self.msg = str(msg)

    def __str__(self):
        if self.msg :
            s = self.msg + ': '
        else :
            s = ''

        return s+self.file

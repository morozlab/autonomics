
class PseudoContig:
    
    def __init__(self, seqID, baseSequence):
        self.leftOverhangLen = 0
        self.rightOverhangLen = 0
        self.positions = []
        self.extensions = []
        self.alignments = {}
        self.id = seqID
        self.baseSeq = baseSequence
        self.baseLen = len(baseSequence)
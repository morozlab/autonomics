
class ContigExtension:
    
    def __init__(self, overhang, direction, ID, sequence, startIndex):
        self.overhang = overhang
        self.direction = direction
        self.extensionID = ID
        self.sequence = sequence
        self.posInDBSeq = startIndex
        
    
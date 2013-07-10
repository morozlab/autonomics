#this class defines a contained for the bases seen at a particular position in a forming contig or scaffold
class BaseWeights:
    
    def __init__(self):
        self.A = 0
        self.T = 0
        self.C = 0
        self.G = 0
        self.N = 0
        self.backbone = 'N'
        self.consensusBase = 'A'
    
        
    
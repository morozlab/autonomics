class Credentials:

    def __init__(self, host, user, passwd):
        self.host = host
        self.user = user
        self.passwd = passwd

    def __getitem__(self, item_name):
        return getattr(self, item_name)
    
    def __setitem__(self, item_name, item_value):
        setattr(self, item_name, item_value)

    def update(self, u, p):
        self.user = u
        self.passwd = p

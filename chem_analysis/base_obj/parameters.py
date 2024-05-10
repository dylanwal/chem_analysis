
class Parameters:

    def __init__(self, **kwargs):
        if kwargs:
            for k, v in kwargs.items():
                setattr(self, k, v)


from chem_analysis.base_obj.parameters import Parameters


class GCParameters(Parameters):
    def __init__(self,
                 *,
                 method: str = None,
                 start_time: int | float = None,
                 end_time: int | float = None,
                 **kwargs
                 ):
        super().__init__(**kwargs)
        self.method = method
        self.start_time = start_time
        self.end_time = end_time

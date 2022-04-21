from functools import wraps


def up_to_date(func):
    @wraps(func)
    def _update(*args, **kwargs):
        if not args[0]._up_to_date:
            args[0]._update()
        return func(*args, **kwargs)
    return _update

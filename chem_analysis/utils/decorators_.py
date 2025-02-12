from functools import wraps
import logging


def try_else_zero(logger_: logging.Logger):
    def _try_else_zero(func):
        @wraps(func)
        def wrap_try_else_zero(*args, **kwargs):
            try:
                return func(*args, **kwargs)
            except Exception as e:
                logger_.warning(f"{func.__name__} had the following error:" + e.args[0])
                return 0
        return wrap_try_else_zero
    return _try_else_zero

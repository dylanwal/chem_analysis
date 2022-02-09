from typing import Union


def sig_figs(number: Union[float, int, None], sig_digit: int = 3) -> Union[int, float, None]:
    """ significant figures
    Given a number return a string rounded to the desired significant digits.
    Parameters
    ----------
    number: float, int
        number you want to reduce significant figures on
    sig_digit: int
        significant digits
    Returns
    -------
    number: int, float
    """
    if isinstance(number, float):
        return float('{:.{p}g}'.format(number, p=sig_digit))
    elif isinstance(number, int):
        return int('{:.{p}g}'.format(number, p=sig_digit))
    elif number is None:
        return None
    else:
        raise TypeError(f"'sig_figs' only accepts int or float. Given: {number} (type: {type(number)}")

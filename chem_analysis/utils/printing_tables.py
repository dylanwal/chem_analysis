from __future__ import annotations
from typing import Sequence
from logging import getLogger

from tabulate import tabulate

from chem_analysis.config import global_config

logger = getLogger(__name__)


def apply_sig_figs(number: float | int | None, sig_digit: int = 3) -> int | float | None:
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
        return float('{:0.{}}'.format(number, sig_digit))
    elif isinstance(number, int):
        return int(float('{:0.{}}'.format(float(number), sig_digit)))
    elif number is None:
        return None
    else:
        raise TypeError(f"'sig_figs' only accepts int or float. Given: {number} (type: {type(number)}")


class StatsTable:

    def __init__(self, rows: list[list], headers: list[str]):
        self.rows = rows
        self.headers = headers

    def __str__(self):
        return self.to_str()

    def __repr__(self):
        return f"rows: {len(self.rows)}, cols: {len(self.headers)}"

    def get(self, header: str) -> list[str]:
        if header not in self.headers:
            from rapidfuzz import process
            top_matches = process.extract(header, self.headers, limit=3)
            top_matches = [m[0] for m in top_matches]
            raise KeyError(f"{type(self).__name__}.get(): '{header}' is not in 'headers'."
                           f"\n\tDid you mean: {' or '.join(top_matches)}")

        index = self.headers.index(header)
        return [row[index] for row in self.rows]

    def join(self, table: StatsTable, include_empty: bool = True):
        headers = self.headers + table.headers
        # remove duplicates while maintaining order
        seen = set()
        self.headers = [x for x in headers if x not in seen and not seen.add(x)]

        if table.rows:
            self.rows += table.rows
        else:
            if include_empty:
                self.rows += [[None] * len(self.headers)]

    def to_str(self,
               limit_to: Sequence[str] = None,
               exclude: Sequence[str] = None,
               sig_figs: int = global_config.sig_fig,
               **kwargs
               ):
        if "tablefmt" not in kwargs:
            kwargs["tablefmt"] = global_config.table_format
        rows = process_rows_to_str(self.rows, sig_figs)

        if limit_to is not None:
            headers, rows = do_limit_to(self.headers, rows, limit_to)
        elif exclude is not None:
            headers, rows = do_exclude(self.headers, rows, limit_to)
        else:
            headers = self.headers

        return tabulate(rows, headers, **kwargs)

    def to_csv_str(self,
                   limit_to: Sequence[str] = None,
                   exclude: Sequence[str] = None,
                   with_headers: bool = True,
                   sig_figs: int = global_config.sig_fig
                   ) -> str:
        rows = process_rows_to_str(self.rows, sig_figs)

        if limit_to is not None:
            headers, rows = do_limit_to(self.headers, rows, limit_to)
        elif exclude is not None:
            headers, rows = do_exclude(self.headers, rows, limit_to)
        else:
            headers = self.headers

        if with_headers:
            inner_strings = [",".join(headers)]
        else:
            inner_strings = []

        inner_strings += [",".join(map(str, row)) for row in rows]
        return "\n".join(inner_strings)

    @classmethod
    def from_dict(cls, dict_: dict) -> StatsTable:
        headers = list(dict_.keys())
        return StatsTable(rows=[values_from_dict(dict_, headers, 0)], headers=["peak"] + headers)

    @classmethod
    def from_list_dicts(cls, list_: list[dict]) -> StatsTable:
        headers = get_headers_from_list_dicts(list_)
        return StatsTable(rows=values_from_list_of_dict(list_, headers), headers=["peak"] + headers)


def process_rows_to_str(rows: list[list], sig_figs: int) -> list[list]:
    rows_ = []
    for row in rows:
        row_ = []
        for v in row:
            row_.append(convert_to_str(v, sig_figs))
        rows_.append(row_)

    return rows_


def convert_to_str(value, sig_figs: int):
    if isinstance(value, float) or isinstance(value, int):
        value = apply_sig_figs(value, sig_figs)

    return value


def get_headers_from_list_dicts(list_) -> list[str]:
    if len(list_) == 0:
        return []

    headers = list(list_[0].keys())
    if len(list_) == 0:
        return headers

    for dict_ in list_[1:]:
        keys = dict_.keys()
        for k in keys:
            if k not in headers:
                headers.append(k)

    return headers


def values_from_dict(dict_: dict, headers: list[str], index: int) -> list:
    values = [index]
    for header in headers:
        values.append(dict_.get(header))

    return values


def values_from_list_of_dict(list_: list[dict], headers: list[str]) -> list:
    rows = []
    for i, dict_ in enumerate(list_):
        rows.append(values_from_dict(dict_, headers, i))

    return rows


def do_limit_to(headers: Sequence[str], rows: list[list[str]], limit_to: Sequence[str]) \
        -> tuple[list[str], list[list[str]]]:
    index = []
    for limit in limit_to:
        if limit not in headers:
            logger.warning(f"'{limit}' is not in headers: {headers}")
            continue
        index.append(headers.index(limit))

    headers_ = [headers[i] for i in index]
    rows_ = []
    for i in range(len(rows)):
        rows_.append([v for ii, v in enumerate(rows[i]) if ii in index])

    return headers_, rows_


def do_exclude(headers: Sequence[str], rows: list[list[str]], exclude: Sequence[str]) \
        -> tuple[list[str], list[list[str]]]:
    rows_ = []
    headers_ = []
    exclude = list(set(exclude))
    for i in range(len(headers)):
        if headers[i] in exclude:
            exclude.remove(headers[i])
            continue
        rows_.append(rows[i])
        headers_.append(headers[i])

    if len(exclude) != 0:
        logger.warning(f"'{exclude}' is not in headers: {headers}")

    return headers_, rows_

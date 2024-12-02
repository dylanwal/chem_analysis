import pathlib
import csv
import logging

from chem_analysis.base_obj.signal_ import Signal
from chem_analysis.analysis.peak import PeakContinuousData
from chem_analysis.analysis.integration.result_integration import ResultIntegration, ResultIntegration2D
from chem_analysis.utils.math import get_slice

logger = logging.getLogger(__name__)


def parse_and_validate_row(row, row_num):
    """Validate a single row."""
    if len(row) != 3:
        raise ValueError(f"Row {row_num} does not have exactly 3 columns: {row}")
    try:
        float(row[1])  # Check if column 2 can be converted to float
        float(row[2])  # Check if column 3 can be converted to float
    except ValueError:
        raise ValueError(f"Row {row_num}, Columns 2 and 3 must be numbers: {row}")


def is_header(row: list[str]) -> bool:
    if len(row) != 3:
        return False
    try:
        float(row[1])
        return False
    except ValueError:
        return True


def load_table(file_path: str | pathlib.Path):
    with open(file_path, mode='r', encoding='utf-8') as file:
        rows = []
        reader = csv.reader(file)
        header = next(reader)
        if not is_header(header):
            parse_and_validate_row(header, 0)
        for i, row in enumerate(csv.reader(file), start=1):
            parse_and_validate_row(row, i)  # Validate each row
            rows.append([row[0], float(row[1]), float(row[2])])

        if len(rows) == 0:
            logger.warning("file is empty.")
            raise ValueError(f"File {file_path} is empty.")

        return rows


def integrate_from_file_trapz(signal: Signal, file: str | pathlib.Path) -> ResultIntegration | ResultIntegration2D:
    # load file
    # format of file must csv: str,float,float  or label,left bound, right bound
    lib = load_table(file)

    result = ResultIntegration(signal=signal)
    for i, row in enumerate(lib):
        lb = row[1]
        ub = row[2]
        slice_ = get_slice(signal.x, lb, ub)
        result.add_peak(
            PeakContinuousData(
                parent=signal,
                slice_=slice_,
                label=row[0],
                id_=i
            )
        )

    if len(result.peaks) == 0:
        logger.warning("No peaks detection during integration.")

    return result

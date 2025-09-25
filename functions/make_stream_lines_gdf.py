""" Function to get supply adjustment specification data from a Supply Adjustment Specifications File"""

import datetime

import pandas as pd

from ..files.get_var import get_var
from ..files.line_to_list import line_to_list
from ..files.skip_until_flag import skip_until_flag

def get_supply_adj_specs(
        supply_adj_specs_file_path: str,
        ncoladj_flag: str = "/ NCOLADJ",
        eob_flag: str="C-----------",
        skip_blocks: int=2
)-> pd.DataFrame:

    """
    Parse supply adjustment specifications from a Supply Adjustment Specifications File.

    This function locates the NCOLADJ variable (number of columns in the supply adjustment specification data file),
    navigates to the start of the supply adjustment specifications table by skipping a configurable number of
    end-of-block markers, and then parses supply adjustment options for each timestep row into a structured table.

    Parameters
    ----------
    supply_adj_specs_file_path : str
        Path to the Supply Adjustment Specifications file (text format).
    ncoladj_flag : str, optional
        Flag used to locate the NCOLADJ value (number of columns in the supply adjustment specifications data file),
        by default "/ NCOLADJ".
    eob_flag : str, optional
        Marker indicating the end of a data block or delimiting headers, by default "C-----------".
    skip_blocks : int, optional
        Number of end-of-block markers to skip after NCOLADJ to reach the supply adjustment specifications table,
        by default 2.

    Returns
    -------
    pandas.DataFrame
        A table with one row per timestep and the following columns:
        - 'date' (datetime): Timestamp.
        - 'kadj_k' (str): Supply adjustment option, for each k in 1 to NCOLADJ.

    Raises
    ------
    ValueError
        If NCOLADJ cannot be read, if the file ends before all timesteps are parsed,
        if a row has fewer than the expected number of fields, if the date format is invalid.
    """

    with open(supply_adj_specs_file_path, "r") as read_file:
        data = read_file.read()

    # List with lines
    lines = data.split("\n")

    # Let's find the number of supply adjustment columns
    try:
        ncoladj, i = get_var(supply_adj_specs_file_path, ncoladj_flag)
        ncoladj = int(ncoladj)

    except Exception as exc:
        raise ValueError(f"Could not read NCOLADJ from flag '{ncoladj_flag}' in {supply_adj_specs_file_path}") from exc

    # Let's skip blocks
    for _ in range(skip_blocks):
        i = skip_until_flag(lines, i, eob_flag)
        i += 1

    # Let's find theindex of the end of the table
    i_eot = skip_until_flag(lines, i, eob_flag)

    records = []
    k = 1
    while i < i_eot:

        line = lines[i]
        line_list = line_to_list(line)

        if len(line_list) < ncoladj + 1:
            raise ValueError(f"supply adjustment row {k} at line {i} has {len(line_list)} fields, "
                             f"expected ≥ {ncoladj + 1}: '{line}'")

        record = {}
        try:
            record["date"] = datetime.datetime.strptime(line_list[0], "%m/%d/%Y_24:00")
        except ValueError:
            raise ValueError(f"Invalid date format in line: {line}")

        for j in range(1, ncoladj + 1):
            record[f"kadj_{j}"] = str(line_list[j])

        records.append(record)

        i += 1
        k += 1

    df = pd.DataFrame(
        records
    )

    return df
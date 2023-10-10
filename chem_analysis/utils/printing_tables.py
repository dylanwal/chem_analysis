
def _length_limit(label: str, limit: int) -> str:
    """ For printing out stats. """
    if len(label) > limit:
        return label[0:limit - 5] + "..."

    return label


def generate_table_from_dict(dict_: dict) -> str:
    text = ""
    # format
    width = int(window / len(headers))
    row_format = ("{:<" + str(width) + "}") * len(headers)

    # headers
    if op_headers:
        headers_ = [_length_limit(head, width) for head in headers.values()]
        text = row_format.format(*headers_) + "\n"
        text = text + "-" * window + "\n"

    # values
    entries = [_length_limit(str(sig_figs(getattr(self, k), 3)), width) for k in headers]
    entries[0] = f"{self.parent.name}: {entries[0]}"  # peak name: peak id
    text = text + row_format.format(*entries) + "\n"

    if op_print:
        print(text)

    return text